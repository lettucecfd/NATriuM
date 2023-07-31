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

ShearLayerStats::ShearLayerStats(CompressibleCFDSolver<3> &solver, std::string outdir, double starting_delta_theta, double starting_Re) :
        DataProcessor<3>(solver),
        m_Re0(starting_Re), m_u(solver.getVelocity()), m_rho(solver.getDensity()),
        m_outDir(outdir), m_filename(scalaroutfile(solver.getConfiguration()->getOutputDirectory())),
        m_vectorfilename(vectoroutfile(solver.getConfiguration()->getOutputDirectory())),
        m_currentDeltaTheta(starting_delta_theta), m_currentDeltaOmega(0.41), m_b11(0), m_b22(0), m_b12(0) {

    m_yCoordsUpToDate = false;
    m_nofCoordinates = 0;

    if (solver.getIterationStart() > 0) {
        if (is_MPI_rank_0()) {
            m_tableFile = boost::make_shared<std::fstream>(m_filename, std::fstream::out | std::fstream::app);
            m_vectorFile = boost::make_shared<std::fstream>(m_vectorfilename, std::fstream::out | std::fstream::app);
        }
    } else {
        if (is_MPI_rank_0()) {
            m_tableFile = boost::make_shared<std::fstream>(m_filename, std::fstream::out);
            m_vectorFile = boost::make_shared<std::fstream>(m_vectorfilename, std::fstream::out);
        }
    }
    if (is_MPI_rank_0()) {
        *m_tableFile << "it ";
        *m_tableFile << "t ";
        *m_tableFile << "deltaTheta ";
        *m_tableFile << "deltaOmega ";
        *m_tableFile << "deltaOmegaDot ";
        *m_tableFile << "m_b11 ";
        *m_tableFile << "m_b22 ";
        *m_tableFile << "m_b12 ";
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
            const std::vector<dealii::Point<3> >& quad_points = fe_values.get_quadrature_points();
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
//    size_t j;
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
    vector<double> ux_average;
    vector<double> uy_average;
    vector<double> uz_average;
    vector<double> rhoux_average;
    vector<double> rhou11_average;
    vector<double> rhou22_average;
    vector<double> rhou33_average;
    vector<double> rhou12_average;
    vector<double> rho_average;
    vector<double> ux_favre;
    vector<double> momentumthickness_integrand;
    vector<double> umag_average;
    vector<double> dUdy_abs;
    // resize to fit length of y
    rhoux_average.resize(m_nofCoordinates);
    rho_average.resize(m_nofCoordinates);
    ux_favre.resize(m_nofCoordinates);
    momentumthickness_integrand.resize(m_nofCoordinates);
    umag_average.resize(m_nofCoordinates);
    dUdy_abs.resize(m_nofCoordinates);
    m_R11.resize(m_nofCoordinates);
    m_R22.resize(m_nofCoordinates);
    m_R33.resize(m_nofCoordinates);
    m_R12.resize(m_nofCoordinates);
    rhou11_average.resize(m_nofCoordinates);
    rhou22_average.resize(m_nofCoordinates);
    rhou33_average.resize(m_nofCoordinates);
    rhou12_average.resize(m_nofCoordinates);
    ux_average.resize(m_nofCoordinates);
    uy_average.resize(m_nofCoordinates);
    uz_average.resize(m_nofCoordinates);

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
                ux_average.at(y_ind) += m_u.at(0)(dof_ind);
                uy_average.at(y_ind) += m_u.at(1)(dof_ind);
                uz_average.at(y_ind) += m_u.at(2)(dof_ind);
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */

    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // add n_points(iy), momentumthickness_integrand(iy), ux_favre(iy), rho*ux(iy), and rho(iy) from different MPI processes
        number.at(iy) = dealii::Utilities::MPI::sum(number.at(iy), MPI_COMM_WORLD);
        // average over number of points at y
        if (number.at(iy) != 0) {
            rhoux_average.at(iy) = dealii::Utilities::MPI::sum(rhoux_average.at(iy), MPI_COMM_WORLD);
            rho_average.at(iy) = dealii::Utilities::MPI::sum(rho_average.at(iy), MPI_COMM_WORLD);
            umag_average.at(iy) = dealii::Utilities::MPI::sum(umag_average.at(iy), MPI_COMM_WORLD);
            rhoux_average.at(iy) /= number.at(iy);
            rho_average.at(iy) /= number.at(iy);
            umag_average.at(iy) /= number.at(iy);
            ux_average.at(iy) /= number.at(iy);
            uy_average.at(iy) /= number.at(iy);
            uz_average.at(iy) /= number.at(iy);
        } else { nonumbers.push_back(iy); }
    }
    // average of neighboring points if there were no points
    for (unsigned long nonumber : nonumbers) {
        size_t iy = nonumber;
        rhoux_average.at(iy) = 0.5 * (rhoux_average.at(iy + 1) + rhoux_average.at(iy - 1));
        rho_average.at(iy) = 0.5 * (rho_average.at(iy + 1) + rho_average.at(iy - 1));
        umag_average.at(iy) = 0;//0.5 * (umag_average.at(iy + 1) + umag_average.at(iy - 1));
        ux_average.at(iy) = 0;
        uy_average.at(iy) = 0;
        uz_average.at(iy) = 0;
    }

    // calculate flux of ux, uy, uz at all points
    double uxflux, uyflux, uzflux;// loop
    for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> > &quad_points = fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                // calculate fluxes
                uxflux = m_u.at(0)(dof_ind) - ux_average.at(y_ind);
                uyflux = m_u.at(1)(dof_ind) - uy_average.at(y_ind);
                uzflux = m_u.at(2)(dof_ind) - uz_average.at(y_ind);
                // fill value vector
                rhou11_average.at(y_ind) += m_rho(dof_ind)*uxflux*uxflux;
                rhou22_average.at(y_ind) += m_rho(dof_ind)*uyflux*uyflux;
                rhou33_average.at(y_ind) += m_rho(dof_ind)*uzflux*uzflux;
                rhou12_average.at(y_ind) += m_rho(dof_ind)*uxflux*uyflux;
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */

    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // average over number of points at y
        if (number.at(iy) != 0) {
            rhou11_average.at(iy) = dealii::Utilities::MPI::sum(rhou11_average.at(iy), MPI_COMM_WORLD);
            rhou22_average.at(iy) = dealii::Utilities::MPI::sum(rhou22_average.at(iy), MPI_COMM_WORLD);
            rhou33_average.at(iy) = dealii::Utilities::MPI::sum(rhou33_average.at(iy), MPI_COMM_WORLD);
            rhou12_average.at(iy) = dealii::Utilities::MPI::sum(rhou12_average.at(iy), MPI_COMM_WORLD);
            rhou11_average.at(iy) /= number.at(iy);
            rhou22_average.at(iy) /= number.at(iy);
            rhou33_average.at(iy) /= number.at(iy);
            rhou12_average.at(iy) /= number.at(iy);
        } else {
            rhou11_average.at(iy) = rhou11_average.at(iy-1);
            rhou22_average.at(iy) = rhou22_average.at(iy-1);
            rhou33_average.at(iy) = rhou33_average.at(iy-1);
            rhou12_average.at(iy) = rhou12_average.at(iy-1);
        }
    }

    // calculate ux_favre, Rij and momentumthickness_integrand
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        ux_favre.at(iy) = rhoux_average.at(iy) / rho_average.at(iy);
        ux_favre.at(iy) = rhoux_average.at(iy) / rho_average.at(iy);
        ux_favre.at(iy) = rhoux_average.at(iy) / rho_average.at(iy);
        momentumthickness_integrand.at(iy) = rho_average.at(iy) * (1 /* dU/2 */ - ux_favre.at(iy) * (1 /* dU/2 */ + ux_favre.at(iy)));
        m_R11.at(iy) = rhou11_average.at(iy) / rho_average.at(iy);
        m_R22.at(iy) = rhou22_average.at(iy) / rho_average.at(iy);
        m_R33.at(iy) = rhou33_average.at(iy) / rho_average.at(iy);
        m_R12.at(iy) = rhou12_average.at(iy) / rho_average.at(iy);
    }

    // calculate dU/dy
    double dy;
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
//        if (iy == 0) { // forward
//            dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
//            dUdy_abs.at(iy) = abs(umag_average.at(iy + 1) - umag_average.at(iy)) / dy;
//        } else if (iy == m_nofCoordinates-1) { // backward
//            dy = m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1);
//            dUdy_abs.at(iy) = abs(umag_average.at(iy) - umag_average.at(iy - 1)) / dy;
//        } else { // other: central
//            dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1);
//            dUdy_abs.at(iy) = abs(umag_average.at(iy - 1) - 2 * umag_average.at(iy) + umag_average.at(iy + 1)) / (dy * dy);
//        }
        if (iy == m_nofCoordinates - 1) { // backward
            dy = m_yCoordinates.at(iy) - m_yCoordinates.at(iy - 1);
            dUdy_abs.at(iy) = abs(umag_average.at(iy) - umag_average.at(iy - 1)) / dy;
        } else { // forward
            dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
            dUdy_abs.at(iy) = abs(umag_average.at(iy + 1) - umag_average.at(iy)) / dy;
        }
    }

    // integrate along y
    double momentumthickness_integral = 0;
    for (size_t iy = 0; iy < m_nofCoordinates - 1; iy++) {
        double window_size;
        if (iy == 0) { // left side: trapezoidal rule
            window_size = abs(m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy));
            momentumthickness_integral += window_size * 0.5 * (momentumthickness_integrand.at(iy) + momentumthickness_integrand.at(iy + 1));
        } else {
            window_size = 0.5*abs(m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1)); // other: simpson rule
            momentumthickness_integral += window_size * (momentumthickness_integrand.at(iy - 1) + 4 * momentumthickness_integrand.at(iy) + momentumthickness_integrand.at(iy + 1)) / 6;
        }
    }

    // calculate momentum thickness
    m_currentDeltaOmega = 2 /*dU*/ / *max_element(std::begin(dUdy_abs), std::end(dUdy_abs));
    m_ReOmega = 1 /*rho0*/ * m_Re0 * 2 /*dU*/ * m_currentDeltaOmega;

    // calculate vorticity thickness
    m_currentDeltaTheta = momentumthickness_integral * (1. / 4. /* rho0 * dU^2 = 1 * 2*2 */);

    // output to console
    if (is_MPI_rank_0()) {
//        cout << "dUdy_abs: ";
//        for (size_t yi = 0; yi < m_nofCoordinates; yi++) {
//            cout << dUdy_abs.at(yi) << ",";
//        } cout << endl;
        cout << "IT: " << m_solver.getIteration()
            << ", t: " << m_solver.getTime()
            << ", delta_Theta: " << m_currentDeltaTheta
            << ", delta_Omega: " << m_currentDeltaOmega
            << ", delta_Omega: " << m_ReOmega
            << ", b11: " << m_b11
            << ", b22: " << m_b22
            << ", b12: " << m_b12
            << endl;
    }
}

void ShearLayerStats::write() {
    if (is_MPI_rank_0()) {
        *m_tableFile << this->m_solver.getIteration() << " "
                     << m_solver.getTime() << " "
                     << m_currentDeltaTheta << " "
                     << m_currentDeltaOmega << " "
                     << m_deltaThetaGrowth << " "
                     << m_b11 << " "
                     << m_b22 << " "
                     << m_b12 << " "
        << endl;

        *m_vectorFile << "IT:" << this->m_solver.getIteration();
        *m_vectorFile << "y: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *m_vectorFile << m_yCoordinates.at(iy) << " ";
        } *m_vectorFile << endl;
        *m_vectorFile << "R11: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *m_vectorFile << m_R11.at(iy) << " ";
        } *m_vectorFile << endl;
        *m_vectorFile << "R22: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *m_vectorFile << m_R22.at(iy) << " ";
        } *m_vectorFile << endl;
        *m_vectorFile << "R33: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *m_vectorFile << m_R33.at(iy) << " ";
        } *m_vectorFile << endl;
        *m_vectorFile << "R12: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *m_vectorFile << m_R12.at(iy) << " ";
        } *m_vectorFile << endl;
    }
}

ShearLayerStats::~ShearLayerStats() = default;

} /* namespace natrium */
