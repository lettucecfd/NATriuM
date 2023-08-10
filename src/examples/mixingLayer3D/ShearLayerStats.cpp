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
    nround = 10000;
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
        *m_tableFile << "deltaThetaDot ";
        *m_tableFile << "deltaOmega ";
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
                y_coords.insert(floor(quad_points.at(i)(1)*nround)/nround);
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
    } // cout << endl;
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
    vector<double> ux_ReAverage(m_nofCoordinates,0), uy_ReAverage(m_nofCoordinates,0),
    uz_ReAverage(m_nofCoordinates,0), umag_average(m_nofCoordinates,0), rho_average(m_nofCoordinates,0),
    rhoux_average(m_nofCoordinates,0), rhouy_average(m_nofCoordinates,0), rhouz_average(m_nofCoordinates,0);
    // initialize vectors for combined sizes
    vector<double> rhou11(m_nofCoordinates,0), rhou22(m_nofCoordinates,0), rhou33(m_nofCoordinates,0), rhou12(m_nofCoordinates,0),
    b11vec(m_nofCoordinates,0), b22vec(m_nofCoordinates,0), b12vec(m_nofCoordinates,0),
    ux_favre(m_nofCoordinates,0), uy_favre(m_nofCoordinates,0), uz_favre(m_nofCoordinates,0),
    momentumthickness_integrand(m_nofCoordinates,0), growthrate_integrand(m_nofCoordinates,0),
    dUdy_abs(m_nofCoordinates,0), uxFavreDy(m_nofCoordinates, 0);
    // resize to fit length of y
    m_R11.resize(m_nofCoordinates);
    m_R22.resize(m_nofCoordinates);
    m_R33.resize(m_nofCoordinates);
    m_R12.resize(m_nofCoordinates);
    m_K.resize(m_nofCoordinates);
    vector<size_t> number(m_nofCoordinates,0);
    vector<size_t> nonumbers(0);

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
                y = floor(quad_points.at(i)(1)*nround)/nround;
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                number.at(y_ind) += 1;
                // fill value vector
                rhoux_average.at(y_ind) += m_rho(dof_ind) * m_u.at(0)(dof_ind); // rho ux
                rhouy_average.at(y_ind) += m_rho(dof_ind) * m_u.at(1)(dof_ind); // rho uy
                rhouz_average.at(y_ind) += m_rho(dof_ind) * m_u.at(2)(dof_ind); // rho uz
                rho_average.at(y_ind) += m_rho(dof_ind);
                umag_average.at(y_ind) += sqrt(pow(m_u.at(0)(dof_ind), 2) + pow(m_u.at(1)(dof_ind), 2) + pow(m_u.at(2)(dof_ind), 2));
                ux_ReAverage.at(y_ind) += m_u.at(0)(dof_ind);
                uy_ReAverage.at(y_ind) += m_u.at(1)(dof_ind);
                uz_ReAverage.at(y_ind) += m_u.at(2)(dof_ind);
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */

    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // add n_points(iy), rho*ux(iy), rho(iy), etc. from different MPI processes
        number.at(iy) = dealii::Utilities::MPI::sum(number.at(iy), MPI_COMM_WORLD);
        // average over number of points at y
        if (number.at(iy) != 0) {
            rhoux_average.at(iy) = dealii::Utilities::MPI::sum(rhoux_average.at(iy), MPI_COMM_WORLD);
            rhouy_average.at(iy) = dealii::Utilities::MPI::sum(rhouy_average.at(iy), MPI_COMM_WORLD);
            rhouz_average.at(iy) = dealii::Utilities::MPI::sum(rhouz_average.at(iy), MPI_COMM_WORLD);
            rho_average.at(iy) = dealii::Utilities::MPI::sum(rho_average.at(iy), MPI_COMM_WORLD);
            umag_average.at(iy) = dealii::Utilities::MPI::sum(umag_average.at(iy), MPI_COMM_WORLD);
            ux_ReAverage.at(iy) = dealii::Utilities::MPI::sum(ux_ReAverage.at(iy), MPI_COMM_WORLD);
            uy_ReAverage.at(iy) = dealii::Utilities::MPI::sum(uy_ReAverage.at(iy), MPI_COMM_WORLD);
            uz_ReAverage.at(iy) = dealii::Utilities::MPI::sum(uz_ReAverage.at(iy), MPI_COMM_WORLD);
            rhoux_average.at(iy) /= number.at(iy);
            rhouy_average.at(iy) /= number.at(iy);
            rhouz_average.at(iy) /= number.at(iy);
            rho_average.at(iy) /= number.at(iy);
            umag_average.at(iy) /= number.at(iy);
            ux_ReAverage.at(iy) /= number.at(iy);
            uy_ReAverage.at(iy) /= number.at(iy);
            uz_ReAverage.at(iy) /= number.at(iy);
        } else { nonumbers.push_back(iy); }
    }
    // average of neighboring points if there were no points
    for (size_t iy : nonumbers) {
        rhoux_average.at(iy) = 0.5 * (rhoux_average.at(iy + 1) + rhoux_average.at(iy - 1));
        rhouy_average.at(iy) = 0.5 * (rhouy_average.at(iy + 1) + rhouy_average.at(iy - 1));
        rhouz_average.at(iy) = 0.5 * (rhouz_average.at(iy + 1) + rhouz_average.at(iy - 1));
        rho_average.at(iy) = 0.5 * (rho_average.at(iy + 1) + rho_average.at(iy - 1));
        umag_average.at(iy) = 0.5 * (umag_average.at(iy + 1) + umag_average.at(iy - 1));
        ux_ReAverage.at(iy) = 0.5 * (ux_ReAverage.at(iy + 1) + ux_ReAverage.at(iy - 1));
        uy_ReAverage.at(iy) = 0.5 * (uy_ReAverage.at(iy + 1) + uy_ReAverage.at(iy - 1));
        uz_ReAverage.at(iy) = 0.5 * (uz_ReAverage.at(iy + 1) + uz_ReAverage.at(iy - 1));
    }

    // calculate flux of ux, uy, uz at all points
    double uxflux_Re, uyflux_Re, uzflux_Re;// loop
    for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> > &quad_points = fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y = floor(quad_points.at(i)(1)*nround)/nround;
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                // calculate fluxes
                uxflux_Re = m_u.at(0)(dof_ind) - ux_ReAverage.at(y_ind);
                uyflux_Re = m_u.at(1)(dof_ind) - uy_ReAverage.at(y_ind);
                uzflux_Re = m_u.at(2)(dof_ind) - uz_ReAverage.at(y_ind);
                // fill value vector
                rhou11.at(y_ind) += m_rho(dof_ind) * uxflux_Re * uxflux_Re;
                rhou22.at(y_ind) += m_rho(dof_ind) * uyflux_Re * uyflux_Re;
                rhou33.at(y_ind) += m_rho(dof_ind) * uzflux_Re * uzflux_Re;
                rhou12.at(y_ind) += m_rho(dof_ind) * uxflux_Re * uyflux_Re;
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */
    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // average over number of points at y
        if (number.at(iy) != 0) {
            rhou11.at(iy) = dealii::Utilities::MPI::sum(rhou11.at(iy), MPI_COMM_WORLD);
            rhou22.at(iy) = dealii::Utilities::MPI::sum(rhou22.at(iy), MPI_COMM_WORLD);
            rhou33.at(iy) = dealii::Utilities::MPI::sum(rhou33.at(iy), MPI_COMM_WORLD);
            rhou12.at(iy) = dealii::Utilities::MPI::sum(rhou12.at(iy), MPI_COMM_WORLD);
            rhou11.at(iy) /= number.at(iy);
            rhou22.at(iy) /= number.at(iy);
            rhou33.at(iy) /= number.at(iy);
            rhou12.at(iy) /= number.at(iy);
//        } else {
//            rhou11.at(iy) = rhou11.at(iy - 1);
//            rhou22.at(iy) = rhou22.at(iy - 1);
//            rhou33.at(iy) = rhou33.at(iy - 1);
//            rhou12.at(iy) = rhou12.at(iy - 1);
        }
    }
    // average of neighboring points if there were no points
    for (size_t iy : nonumbers) {
        rhou11.at(iy) = 0.5 * (rhou11.at(iy + 1) + rhou11.at(iy - 1));
        rhou22.at(iy) = 0.5 * (rhou22.at(iy + 1) + rhou22.at(iy - 1));
        rhou33.at(iy) = 0.5 * (rhou33.at(iy + 1) + rhou33.at(iy - 1));
        rhou12.at(iy) = 0.5 * (rhou12.at(iy + 1) + rhou12.at(iy - 1));
    }

    if (is_MPI_rank_0()) {
        // calculate ui_favre, Rij and momentumthickness_integrand
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            ux_favre.at(iy) = rhoux_average.at(iy) / rho_average.at(iy);
            uy_favre.at(iy) = rhouy_average.at(iy) / rho_average.at(iy);
            uz_favre.at(iy) = rhouz_average.at(iy) / rho_average.at(iy);
            m_R11.at(iy) = rhou11.at(iy) / rho_average.at(iy);
            m_R22.at(iy) = rhou22.at(iy) / rho_average.at(iy);
            m_R33.at(iy) = rhou33.at(iy) / rho_average.at(iy);
            m_R12.at(iy) = rhou12.at(iy) / rho_average.at(iy);
            momentumthickness_integrand.at(iy) = rho_average.at(iy) * (1 /* dU/2 */ - ux_favre.at(iy) * (1 /* dU/2 */ + ux_favre.at(iy)));
            // calculate turbulent kinetic energy
            m_K.at(iy) = (m_R11.at(iy) + m_R22.at(iy) + m_R33.at(iy))/2; // https://www.osti.gov/pages/servlets/purl/1580489
            // calculate anisotropy tensor along y
            b11vec.at(iy) = (m_R11.at(iy) - 2./3*m_K.at(iy)*1)/(2*m_K.at(iy));
            b22vec.at(iy) = (m_R22.at(iy) - 2./3*m_K.at(iy)*1)/(2*m_K.at(iy));
            b12vec.at(iy) = (m_R12.at(iy) - 2./3*m_K.at(iy)*0)/(2*m_K.at(iy));
        }

        // calculate dU/dy
        double dy;
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            if (iy == 0) { // forward
                dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
                dUdy_abs.at(iy) = (umag_average.at(iy + 1) - umag_average.at(iy)) / dy;
                uxFavreDy.at(iy) = (ux_favre.at(iy + 1) - ux_favre.at(iy)) / dy;
            } else if (iy == m_nofCoordinates-1) { // backward
                dy = m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1);
                dUdy_abs.at(iy) = (umag_average.at(iy) - umag_average.at(iy - 1)) / dy;
                uxFavreDy.at(iy) = (ux_favre.at(iy) - ux_favre.at(iy - 1)) / dy;
            } else { // other: central
                dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1);
                dUdy_abs.at(iy) = (umag_average.at(iy + 1) -  umag_average.at(iy - 1)) / dy;
                uxFavreDy.at(iy) = (ux_favre.at(iy + 1) - ux_favre.at(iy - 1)) / dy;
            }
            dUdy_abs.at(iy) = abs(dUdy_abs.at(iy));
            growthrate_integrand.at(iy) = rho_average.at(iy) * m_R12.at(iy) * uxFavreDy.at(iy);
        }

        // integrate along y
        double momentumthickness_integral = 0, growthrate_integral = 0;
        double window_size;
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            if (iy == 0) { // left side: trapezoidal rule
                window_size = (m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy));
                momentumthickness_integral += window_size * 0.5 * (momentumthickness_integrand.at(iy) + momentumthickness_integrand.at(iy + 1));
//                cout << momentumthickness_integral << " ";
                growthrate_integral += window_size * 0.5 * (growthrate_integrand.at(iy) + growthrate_integrand.at(iy + 1));
//                cout << growthrate_integral << endl;
            } else if (iy == m_nofCoordinates - 1) {
                window_size = (m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1));
                momentumthickness_integral += window_size * 0.5 * (momentumthickness_integrand.at(iy) + momentumthickness_integrand.at(iy -1));
                growthrate_integral += window_size * 0.5 * (growthrate_integrand.at(iy) + growthrate_integrand.at(iy - 1));
            } else { // other: simpson rule
                window_size = 0.5*(m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1));
                momentumthickness_integral += window_size * (momentumthickness_integrand.at(iy - 1) + 4 * momentumthickness_integrand.at(iy) + momentumthickness_integrand.at(iy + 1)) / 6;
                growthrate_integral += window_size * (growthrate_integrand.at(iy - 1) + 4 * growthrate_integrand.at(iy) + growthrate_integrand.at(iy + 1)) / 6;
            }
        }

        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            // calculate anisotropy tensor elements
            y = m_yCoordinates.at(iy);
            if (y > -0.093 && y < 0.093) { // simpson rule
                window_size = 0.5 * abs(m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1));
                m_b11 += window_size * (b11vec.at(iy - 1) + 4 * b11vec.at(iy) + b11vec.at(iy + 1)) / 6;
                m_b22 += window_size * (b22vec.at(iy - 1) + 4 * b22vec.at(iy) + b22vec.at(iy + 1)) / 6;
                m_b12 += window_size * (b12vec.at(iy - 1) + 4 * b12vec.at(iy) + b12vec.at(iy + 1)) / 6;
            }
        }

        // calculate vorticity thickness
        m_currentDeltaOmega = 2 /*dU*/ / *max_element(std::begin(dUdy_abs), std::end(dUdy_abs));
        m_ReOmega = 1 /*rho0*/ * 2 /*dU*/ * m_Re0 * m_currentDeltaOmega;

        // calculate momentum thickness
        m_currentDeltaTheta = momentumthickness_integral * (1. / 4. /* rho0 * dU^2 = 1 * 2*2 */);
        m_deltaThetaGrowth = growthrate_integral * (-2) / (8 /* rh0*dU^3 = 1 * 2*2*2 */);

        // output to console
//        cout << "dUdy_abs: ";
//        for (size_t yi = 0; yi < m_nofCoordinates; yi++) {
//            cout << dUdy_abs.at(yi) << ",";
//        } cout << endl;
        cout << "IT: " << m_solver.getIteration()
            << ", t: " << m_solver.getTime()
            << ", delta_Theta: " << m_currentDeltaTheta
            << ", growth_rate: " << m_deltaThetaGrowth
            << ", delta_Omega: " << m_currentDeltaOmega
            << ", Re_Omega: " << m_ReOmega
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
                     << m_deltaThetaGrowth << " "
                     << m_currentDeltaOmega << " "
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
