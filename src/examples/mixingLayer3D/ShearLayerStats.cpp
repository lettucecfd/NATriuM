/*
 * ShearLayerStats.cpp
 *
 *  Created on: 06.10.2022
 *      Author: dwilde3m
 */

#include "ShearLayerStats.h"

#include "mpi.h"
#include <utility>
#include "deal.II/base/mpi.h"
//#include "natrium/collision_advanced/AuxiliaryCollisionFunctions.h"

namespace natrium {

ShearLayerStats::ShearLayerStats(CompressibleCFDSolver<3> &solver, std::string outdir, double starting_delta_theta) :
DataProcessor<3>(solver), m_u(solver.getVelocity()), m_rho(solver.getDensity()),
m_outDir(outdir), m_filename(outfile(solver.getConfiguration()->getOutputDirectory())),
m_currentRho(1.0), m_currentDeltaTheta(starting_delta_theta), m_currentRhoUx(1.0), m_currentUxFavre(1.0) {

//    m_DeltaTheta.push_back(m_currentDeltaTheta);

    m_names.emplace_back("deltaTheta");
    m_names.emplace_back("Rho");
    m_names.emplace_back("RhoUx");
    m_names.emplace_back("UxFavre");
    m_names.emplace_back("T");
    m_names.emplace_back("dT/dx");
    m_names.emplace_back("dT/dy");
    m_names.emplace_back("dT/dz");
    m_names.emplace_back("Ma_local/Ma_wall");
    m_names.emplace_back("rho*ux");
    m_names.emplace_back("rho*uy");
    m_names.emplace_back("rho*uz");

    m_yCoordsUpToDate = false;
    m_nofCoordinates = 0;

    if (solver.getIterationStart() > 0) {
        if (is_MPI_rank_0()) {
            m_tableFile = boost::make_shared<std::fstream>(m_filename, std::fstream::out | std::fstream::app);
        }
    } else {
        if (is_MPI_rank_0()) {
            m_tableFile = boost::make_shared<std::fstream>(m_filename, std::fstream::out);
        }
    }

    *m_tableFile << "it ";
    *m_tableFile << "t    ";
    *m_tableFile << "deltaTheta ";
    *m_tableFile << "rho    ";
    *m_tableFile << "rhoUx  ";
    *m_tableFile << "UxFavre ";
    *m_tableFile << endl;
}

bool ShearLayerStats::isMYCoordsUpToDate() const {
    return m_yCoordsUpToDate;
}

void ShearLayerStats::updateYValues() {
    boost::shared_ptr<AdvectionOperator<3> > advection =
            m_solver.getAdvectionOperator();
    m_yCoordsUpToDate = true;

    //////////////////////////
    // Calculate y values ////
    //////////////////////////
    std::set<double, own_double_less> y_coords;

    const dealii::UpdateFlags update_flags = dealii::update_quadrature_points;
    const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
    dealii::FEValues<3> fe_values(advection->getMapping(),
                                  *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
    // loop
    typename dealii::DoFHandler<3>::active_cell_iterator cell =
            dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            // get y coordinates
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> >& quad_points =
                    fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y_coords.insert(quad_points.at(i)(1));
            }
        }
    }

    //communicate list of y values
    size_t n_ycoords = y_coords.size();
    size_t max_nycoords = dealii::Utilities::MPI::min_max_avg(n_ycoords,
                                                              MPI_COMM_WORLD).max;
    auto *sendbuf = (double*) malloc(max_nycoords * sizeof(double));

    // fill send buffer with y coordinates
    size_t i = 0;
    size_t j;
    for (double y_coord : y_coords) {
        sendbuf[i] = y_coord;
        i++;
    }
    for (; i < max_nycoords; i++) {
        sendbuf[i] = *y_coords.begin();
    }
    // allocate read buffer
    double *recvbuf = (double*) malloc(dealii::Utilities::MPI::n_mpi_processes(
            MPI_COMM_WORLD) * max_nycoords * sizeof(double));
    // allgather: distribute all ycoords to all processors
    MPI_Allgather(sendbuf, max_nycoords,
                  MPI_DOUBLE, recvbuf, max_nycoords,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // remove double entries
    std::set<double> y_coords_gathered;
    for (i = 0;
         i < max_nycoords * dealii::Utilities::MPI::n_mpi_processes(
                 MPI_COMM_WORLD); i++) {
        y_coords_gathered.insert(recvbuf[i]);
    }
    // fill member variables
    i = 0;
    for (double it : y_coords_gathered) {
        m_yCoordinates.push_back(it);
        m_yCoordinateToIndex.insert(std::make_pair(it, i));
        i++;
    }
    m_nofCoordinates = i;
    // free
    free(sendbuf);
    free(recvbuf);
}

void ShearLayerStats::apply() {
    if (!isMYCoordsUpToDate()) {
        updateYValues();
    }
	if (m_solver.getIteration() % m_solver.getConfiguration()->getOutputShearLayerInterval() == 0) {
        calculateRhoU();
        rescaleDensity();
        write();
	}
}

void ShearLayerStats::calculateRhoU() {
    //////////////////////////
    // Calculate averages ////
    //////////////////////////
    // initialize vectors for averages along x and z
    vector<double> rhoux_average;
    vector<double> rho_average;
    vector<double> ux_favre;
    vector<double> integrand;
    // resize to fit length of y
    rhoux_average.resize(m_nofCoordinates);
    rho_average.resize(m_nofCoordinates);
    ux_favre.resize(m_nofCoordinates);
    integrand.resize(m_nofCoordinates);

    vector<size_t> number;
    number.resize(m_nofCoordinates);

    // don't know what I do here, but it worked for turbulent channel
    boost::shared_ptr<AdvectionOperator<3> > advection = m_solver.getAdvectionOperator();
    const dealii::UpdateFlags update_flags = dealii::update_quadrature_points
                                             | dealii::update_gradients;
    const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
    dealii::FEValues<3> fe_values(advection->getMapping(),
                                  *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
    size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
    std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);

    // loop
    typename dealii::DoFHandler<3>::active_cell_iterator cell =
            dof_handler.begin_active(), endc = dof_handler.end();
    double y;
    size_t y_ind;
    size_t dof_ind;
    // get rho and rhoux Reynolds averaged - i.e. along x and z
    for (; cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> >& quad_points =
                    fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                number.at(y_ind) += 1;
                // fill value vector
                // add to averages:
                rhoux_average.at(y_ind) += m_rho(dof_ind) * m_u.at(0)(dof_ind);				// rho u
                rho_average.at(y_ind) += m_rho(dof_ind);
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */

    // communicate
    for (size_t i = 0; i < m_nofCoordinates; i++) {
        number.at(i) = dealii::Utilities::MPI::sum(number.at(i), MPI_COMM_WORLD);
        ux_favre.at(i) = rhoux_average.at(i) / rho_average.at(i);
        integrand.at(i) = rho_average.at(i) * (1 /* dU / 2 */ - ux_favre.at(i) * (1 /* dU / 2 */ + ux_favre.at(i)));
        ux_favre.at(i) = dealii::Utilities::MPI::sum(ux_favre.at(i), MPI_COMM_WORLD);
        rhoux_average.at(i) = dealii::Utilities::MPI::sum(rhoux_average.at(i), MPI_COMM_WORLD);
        rho_average.at(i) = dealii::Utilities::MPI::sum(rho_average.at(i), MPI_COMM_WORLD);
        ux_favre.at(i) /= number.at(i);
        rhoux_average.at(i) /= number.at(i);
        rho_average.at(i) /= number.at(i);
    }
    // integrate along y
    double integral = 0;
    double rho_int = 0;
    double rhoux_int = 0;
    double ux_favre_int = 0;
    double interval_length = 0;
    for (size_t i = 0; i < m_nofCoordinates-1; i++) {
        double window_size = std::abs( m_yCoordinates.at(i+1) -m_yCoordinates.at(i));
        integral += window_size*0.5*(integrand.at(i)+integrand.at(i+1));
        rho_int += window_size*0.5*(rho_average.at(i)+rho_average.at(i+1));
        rhoux_int += window_size*0.5*(rhoux_average.at(i)+rhoux_average.at(i+1));
        ux_favre_int += window_size*0.5*(ux_favre.at(i)+ux_favre.at(i+1));
        interval_length += window_size;
    }
    m_currentRho = rho_int / interval_length;
    m_currentRhoUx = rhoux_int / interval_length;
    m_currentUxFavre = ux_favre_int / interval_length;
//    double lastDeltaTheta = m_DeltaTheta.back();
//    double lastTime = m_Time.back();
    double currentTime = m_solver.getTime();
//    m_Time.push_back(currentTime);
    m_currentDeltaTheta = integral * (1. / 4. /* rho0 * dU^2 = 1 * 2*2 */);
//    m_DeltaTheta_diff = m_currentDeltaTheta - lastDeltaTheta;
//    m_DeltaTheta.push_back(m_currentDeltaTheta);
}

void ShearLayerStats::write() {
    if (is_MPI_rank_0()) {
        *m_tableFile << this->m_solver.getIteration() << " ";
        *m_tableFile << this->m_solver.getTime() << " ";
        *m_tableFile << m_currentDeltaTheta << " ";
        *m_tableFile << m_currentRho << " ";
        *m_tableFile << m_currentRhoUx << " ";
        *m_tableFile << m_currentUxFavre << " ";
        *m_tableFile << endl;
    }
}

void ShearLayerStats::rescaleDensity() {
    m_solver.scaleF(1.0/m_currentRho);
}

ShearLayerStats::~ShearLayerStats() = default;

} /* namespace natrium */
