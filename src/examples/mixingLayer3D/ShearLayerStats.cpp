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
m_currentRho(1.0), m_currentDeltaTheta(starting_delta_theta), m_currentDeltaOmega(0.41), m_currentRhoUx(1.0), m_currentUxFavre(1.0),
m_currentTime(0.0) {

//    m_DeltaTheta.push_back(m_currentDeltaTheta);

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
    if (is_MPI_rank_0()) {
        *m_tableFile << "it ";
        *m_tableFile << "t     ";
        *m_tableFile << "deltaTheta ";
        *m_tableFile << "deltaOmega ";
        *m_tableFile << "t*dU/DT0  ";
        *m_tableFile << "DT/DT0 ";
//    *m_tableFile << "rho    ";
//    *m_tableFile << "rhoUx  ";
//    *m_tableFile << "UxFavre ";
        *m_tableFile << endl;
    }
}

bool ShearLayerStats::isMYCoordsUpToDate() const {
    return m_yCoordsUpToDate;
}

void ShearLayerStats::updateYValues() {
    boost::shared_ptr<AdvectionOperator<3> > advection = m_solver.getAdvectionOperator();

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
    // finished
    m_yCoordsUpToDate = true;
}

void ShearLayerStats::apply() {
    if (!isMYCoordsUpToDate()) {
        updateYValues();
    }
	if (m_solver.getIteration() % m_solver.getConfiguration()->getOutputShearLayerInterval() == 0) {
        calculateRhoU();
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
    vector<double> umag_average;
    vector<double> dUdy_abs;
    // resize to fit length of y
    rhoux_average.resize(m_nofCoordinates);
    rho_average.resize(m_nofCoordinates);
    ux_favre.resize(m_nofCoordinates);
    integrand.resize(m_nofCoordinates);
    umag_average.resize(m_nofCoordinates);
    dUdy_abs.resize(m_nofCoordinates);

    vector<size_t> number;
    number.resize(m_nofCoordinates);
    vector<size_t> nonumbers;

    // don't know what I do here, but it worked for turbulent channel
    boost::shared_ptr<AdvectionOperator<3> > advection = m_solver.getAdvectionOperator();
    const dealii::UpdateFlags update_flags = dealii::update_quadrature_points | dealii::update_gradients;
    const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
    dealii::FEValues<3> fe_values(advection->getMapping(),
                                  *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
    size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
    vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);

    // loop
    typename dealii::DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    double y;
    size_t y_ind;
    size_t dof_ind;
    // get rho and rhoux Reynolds averaged - i.e. along x and z
    for (; cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> >& quad_points = fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                number.at(y_ind) += 1;
                // fill value vector
                rhoux_average.at(y_ind) += m_rho(dof_ind) * m_u.at(0)(dof_ind); // rho ux
                rho_average.at(y_ind) += m_rho(dof_ind);
                umag_average.at(y_ind) += sqrt(pow(m_u.at(0)(dof_ind), 2) + pow(m_u.at(1)(dof_ind), 2) + pow(m_u.at(2)(dof_ind), 2));
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */

    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // add n_points(iy), integrand(iy), ux_favre(iy), rho*ux(iy), and rho(iy) from different MPI processes
        number.at(iy) = dealii::Utilities::MPI::sum(number.at(iy), MPI_COMM_WORLD);
        // average over number of points at y
        if (number.at(iy) != 0) {
            rhoux_average.at(iy) = dealii::Utilities::MPI::sum(rhoux_average.at(iy), MPI_COMM_WORLD);
            rho_average.at(iy) = dealii::Utilities::MPI::sum(rho_average.at(iy), MPI_COMM_WORLD);
            umag_average.at(iy) = dealii::Utilities::MPI::sum(umag_average.at(iy), MPI_COMM_WORLD);
            rhoux_average.at(iy) /= number.at(iy);
            rho_average.at(iy) /= number.at(iy);
            umag_average.at(iy) /= number.at(iy);
        } else { nonumbers.push_back(iy); }
    }
    // average of neighboring points if there were no points
    size_t iy;
    for (size_t i = 0; i < nonumbers.size(); i++) {
        iy = nonumbers.at(i);
        rhoux_average.at(iy) = 0;//0.5 * (rhoux_average.at(iy + 1) + rhoux_average.at(iy - 1));
        rho_average.at(iy) = 1;//0.5 * (rho_average.at(iy + 1) + rho_average.at(iy - 1));
        umag_average.at(iy) = 0;//0.5 * (umag_average.at(iy + 1) + umag_average.at(iy - 1));
    }

    // calculate ux_favre and integrand
    for (size_t yi = 0; yi < m_nofCoordinates; yi++) {
        ux_favre.at(yi) = rhoux_average.at(yi) / rho_average.at(yi);
        integrand.at(yi) = rho_average.at(yi) * (1 /* dU 2 */ - ux_favre.at(yi) * (1 /* dU/2 */ + ux_favre.at(yi)));
//        // ignore row with 0 nodes
//        if (number.at(yi) == 0) {
//            integrand.at(yi) = 0;
//        }
    }
    // calculate dU/dy
    double dy;
    for (size_t yi = 0; yi < m_nofCoordinates; yi++) {
        if (yi == 0) { // left hand
            dy = m_yCoordinates.at(yi + 1) - m_yCoordinates.at(yi);
            dUdy_abs.at(yi) = abs(umag_average.at(yi + 1) - umag_average.at(yi)) / dy;
        } else if (yi == m_nofCoordinates-1) { // right hand
            dy = m_yCoordinates.at(yi) - m_yCoordinates.at(yi-1);
            dUdy_abs.at(yi) = abs(umag_average.at(yi) - umag_average.at(yi - 1)) / dy;
        } else { // other: central
            dy = m_yCoordinates.at(yi + 1) - m_yCoordinates.at(yi - 1);
            if (dy < 1e-10) {
                dUdy_abs.at(yi) = dUdy_abs.at(yi-1);
            } else {
                dUdy_abs.at(yi) = abs(umag_average.at(yi - 1) - 2 * umag_average.at(yi) + umag_average.at(yi + 1)) / (dy * dy);
            }
        }
//        // ignore row with 0 nodes
//        if (number.at(yi) == 0) {
//            dUdy_abs.at(yi) = 0;
//        }
    }

    // integrate along y
    double integral = 0;
    double rho_avg = 0;
    double rhoux_avg = 0;
    double ux_favre_avg = 0;
    double interval_length = 0;
    for (size_t yi = 0; yi < m_nofCoordinates-1; yi++) {
        double window_size;
        if (yi == 0) { // left side: trapezoidal rule
            window_size = abs(m_yCoordinates.at(yi + 1) - m_yCoordinates.at(yi));
            integral += window_size * 0.5 * (integrand.at(yi) + integrand.at(yi + 1));
        } else {
            window_size = 0.5*abs(m_yCoordinates.at(yi + 1) - m_yCoordinates.at(yi-1)); // other: simpson rule
            integral += window_size * (integrand.at(yi - 1) + 4 * integrand.at(yi) + integrand.at(yi + 1)) / 6;
        }
        rho_avg += rho_average.at(yi);
        rhoux_avg += rhoux_average.at(yi);
        ux_favre_avg += ux_favre.at(yi);
        interval_length += window_size;
    }
    rho_avg /= m_nofCoordinates-1;
    rhoux_avg /= m_nofCoordinates-1;
    ux_favre_avg /= m_nofCoordinates-1;
    // calculate vorticity thickness
    m_currentDeltaOmega = 2 /*dU*/ / *max_element(std::begin(dUdy_abs), std::end(dUdy_abs));;

    m_currentRho = rho_avg;
    m_currentRhoUx = rhoux_avg;
    m_currentUxFavre = ux_favre_avg;
    // update deltaTheta
    m_lastDeltaTheta = m_currentDeltaTheta;
    m_currentDeltaTheta = integral * (1. / 4. /* rho0 * dU^2 = 1 * 2*2 */);
    // update time
    m_lastTime = m_currentTime;
    m_currentTime = m_solver.getTime();
    double dt = m_currentTime - m_lastTime;
    // calculate difference
    m_DeltaTheta_diff = (m_currentDeltaTheta - m_lastDeltaTheta) / (dt * (2 /*dU*/ / 0.093 /*DT0*/));
    if (is_MPI_rank_0()) {
//        cout << "dUdy_abs: ";
//        for (size_t yi = 0; yi < m_nofCoordinates; yi++) {
//            cout << dUdy_abs.at(yi) << ",";
//        } cout << endl;
        cout << "IT: " << m_solver.getIteration()
            << ", t: " << m_currentTime
            << ", delta_Theta: " << m_currentDeltaTheta
            << ", delta_Omega: " << m_currentDeltaOmega
            << ", t*dU/DT0: " << m_currentTime*2/0.093
            << ", delta_Theta/DT0: " << m_currentDeltaTheta/0.093
            << endl;
    }
}

void ShearLayerStats::write() {
    if (is_MPI_rank_0()) {
        *m_tableFile << this->m_solver.getIteration() << " ";
        *m_tableFile << m_currentTime << " ";
        *m_tableFile << m_currentDeltaTheta << " ";
        *m_tableFile << m_currentDeltaOmega << " ";
        *m_tableFile << m_currentTime*2/*dU*//0.093 << " ";
        *m_tableFile << m_currentDeltaTheta/0.093 << " ";
//        *m_tableFile << m_DeltaTheta_diff << " ";
//        *m_tableFile << m_lastTime << " ";
//        *m_tableFile << m_lastDeltaTheta << " ";
        *m_tableFile << endl;
    }
}

ShearLayerStats::~ShearLayerStats() = default;

} /* namespace natrium */
