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
        m_initializationfilename(initializationoutfile(solver.getConfiguration()->getOutputDirectory())),
        m_currentDeltaTheta(starting_delta_theta), m_currentDeltaOmega(0.41), m_b11(0), m_b22(0), m_b12(0) {
    nround = pow(10,12); // round coordinates to this magnitude
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
    updateYValues();
    calculateRhoU();
    write();

    // run checks for initialization
    if (is_MPI_rank_0()) {
        m_initializationFile = boost::make_shared<std::fstream>(m_initializationfilename, std::fstream::out | std::fstream::app);
        *m_initializationFile << "y: ";
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            *m_initializationFile << m_yCoordinates.at(iy) << ", ";
        } *m_initializationFile << endl;
        *m_initializationFile << "ux_favre: ";
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            *m_initializationFile << ux_favre.at(iy) << ", ";
        } *m_initializationFile << endl;
        *m_initializationFile << "ux_ReAvg: ";
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            *m_initializationFile << ux_Re.at(iy) << ", ";
        } *m_initializationFile << endl;
        *m_initializationFile << "rho_ReAvg: ";
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            *m_initializationFile << rho_average.at(iy) << ", ";
        } *m_initializationFile << endl;
        *m_initializationFile << "momentumthickness_integrand: ";
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            *m_initializationFile << momentumthickness_integrand.at(iy) << ", ";
        } *m_initializationFile << endl;
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

    // get y coordinates
    const dealii::UpdateFlags update_flags = dealii::update_quadrature_points;
    const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
    dealii::FEValues<3> fe_values(advection->getMapping(), *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
    typename dealii::DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> >& quad_points = fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y_coords.insert(floor(quad_points.at(i)(1)*nround)/nround);
//                y_coords.insert(quad_points.at(i)(1));
            }
        }
    }

    //communicate list of y values
    size_t n_ycoords = y_coords.size();
    size_t max_nycoords = dealii::Utilities::MPI::min_max_avg(n_ycoords, MPI_COMM_WORLD).max;
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
    double *recvbuf = (double*) malloc(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) * max_nycoords * sizeof(double));
    // allgather: distribute all ycoords to all processors
    MPI_Allgather(sendbuf, max_nycoords, MPI_DOUBLE, recvbuf, max_nycoords, MPI_DOUBLE, MPI_COMM_WORLD);

    // remove double entries
    std::set<double> y_coords_gathered;
    for (i = 0;
         i < max_nycoords * dealii::Utilities::MPI::n_mpi_processes( MPI_COMM_WORLD); i++) {
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
//    vector<double> ux_Re(m_nofCoordinates, 0), uy_Re(m_nofCoordinates, 0), uz_Re(m_nofCoordinates, 0);
//    vector<double> umag_average(m_nofCoordinates, 0), rho_average(m_nofCoordinates, 0);
    vector<double> rhoux_average(m_nofCoordinates,0), rhouy_average(m_nofCoordinates,0), rhouz_average(m_nofCoordinates,0);
    // initialize vectors for combined sizes
    vector<double> rhou11(m_nofCoordinates,0), rhou22(m_nofCoordinates,0), rhou33(m_nofCoordinates,0), rhou12(m_nofCoordinates,0),
    b11vec(m_nofCoordinates,0), b22vec(m_nofCoordinates,0), b12vec(m_nofCoordinates,0),
//    ux_favre(m_nofCoordinates,0), uy_favre(m_nofCoordinates,0), uz_favre(m_nofCoordinates,0),
//    momentumthickness_integrand(m_nofCoordinates,0),
    growthrate_integrand(m_nofCoordinates,0);
    // resize to fit length of y
    ux_Re.resize(m_nofCoordinates); uy_Re.resize(m_nofCoordinates); uz_Re.resize(m_nofCoordinates);
    umag_average.resize(m_nofCoordinates); rho_average.resize(m_nofCoordinates);
    ux_favre.resize(m_nofCoordinates); uy_favre.resize(m_nofCoordinates); uz_favre.resize(m_nofCoordinates);
    m_R11.resize(m_nofCoordinates); m_R22.resize(m_nofCoordinates); m_R33.resize(m_nofCoordinates);
    m_R12.resize(m_nofCoordinates);
    m_K.resize(m_nofCoordinates);
    momentumthickness_integrand.resize(m_nofCoordinates);
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
//                y = quad_points.at(i)(1);
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
                ux_Re.at(y_ind) += m_u.at(0)(dof_ind);
                uy_Re.at(y_ind) += m_u.at(1)(dof_ind);
                uz_Re.at(y_ind) += m_u.at(2)(dof_ind);
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
            ux_Re.at(iy) = dealii::Utilities::MPI::sum(ux_Re.at(iy), MPI_COMM_WORLD);
            uy_Re.at(iy) = dealii::Utilities::MPI::sum(uy_Re.at(iy), MPI_COMM_WORLD);
            uz_Re.at(iy) = dealii::Utilities::MPI::sum(uz_Re.at(iy), MPI_COMM_WORLD);
            rhoux_average.at(iy) /= number.at(iy);
            rhouy_average.at(iy) /= number.at(iy);
            rhouz_average.at(iy) /= number.at(iy);
            rho_average.at(iy) /= number.at(iy);
            umag_average.at(iy) /= number.at(iy);
            ux_Re.at(iy) /= number.at(iy);
            uy_Re.at(iy) /= number.at(iy);
            uz_Re.at(iy) /= number.at(iy);
        } else nonumbers.push_back(iy);
    }
    // average of neighboring points if there were no points
    int iy_u, iy_l;
    float dy_u, dy_l;
    float window_size;
    for (auto iy : nonumbers) {
        if (iy == 0) {iy_l = 0; iy_u = iy + 1;}
        else if (iy == m_nofCoordinates - 1) {iy_l = iy - 1; iy_u = iy;}
        else {iy_l = iy - 1; iy_u = iy + 1;}
        dy_u = m_yCoordinates.at(iy_u) - m_yCoordinates.at(iy);
        dy_l = m_yCoordinates.at(iy) - m_yCoordinates.at(iy_l);
        window_size = dy_u + dy_l;
        rhoux_average.at(iy) = 1 / window_size * (rhoux_average.at(iy_u) * dy_u + rhoux_average.at(iy_l) * dy_l);
        rhouy_average.at(iy) = 1 / window_size * (rhouy_average.at(iy_u) * dy_u + rhouy_average.at(iy_l) * dy_l);
        rhouz_average.at(iy) = 1 / window_size * (rhouz_average.at(iy_u) * dy_u + rhouz_average.at(iy_l) * dy_l);
        rho_average.at(iy) = 1 / window_size * (rho_average.at(iy_u) * dy_u + rho_average.at(iy_l) * dy_l);
        umag_average.at(iy) = 1 / window_size * (umag_average.at(iy_u) * dy_u + umag_average.at(iy_l) * dy_l);
        ux_Re.at(iy) = 1 / window_size * (ux_Re.at(iy_u) * dy_u + ux_Re.at(iy_l) * dy_l);
        uy_Re.at(iy) = 1 / window_size * (uy_Re.at(iy_u) * dy_u + uy_Re.at(iy_l) * dy_l);
        uz_Re.at(iy) = 1 / window_size * (uz_Re.at(iy_u) * dy_u + uz_Re.at(iy_l) * dy_l);
    }
    auto [minU, maxU] = std::minmax_element(begin(ux_Re), end(ux_Re));
    auto dU_calculated = *maxU - *minU;

    // ui_favre and momentumthickness_integrand
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        ux_favre.at(iy) = rhoux_average.at(iy) / rho_average.at(iy);
        uy_favre.at(iy) = rhouy_average.at(iy) / rho_average.at(iy);
        uz_favre.at(iy) = rhouz_average.at(iy) / rho_average.at(iy);
        momentumthickness_integrand.at(iy) = rho_average.at(iy) * (dU_calculated / 2 - ux_favre.at(iy)) * (dU_calculated / 2 + ux_favre.at(iy));
    }

    // calculate flux of ux, uy, uz at all points
    double uxflux_Fa, uyflux_Fa, uzflux_Fa;// loop
    for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> > &quad_points = fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y = floor(quad_points.at(i)(1)*nround)/nround;
//                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                // calculate fluxes
                uxflux_Fa = m_u.at(0)(dof_ind) - ux_favre.at(y_ind);
                uyflux_Fa = m_u.at(1)(dof_ind) - uy_favre.at(y_ind);
                uzflux_Fa = m_u.at(2)(dof_ind) - uz_favre.at(y_ind);
                // fill value vector
                rhou11.at(y_ind) += m_rho(dof_ind) * uxflux_Fa * uxflux_Fa;
                rhou22.at(y_ind) += m_rho(dof_ind) * uyflux_Fa * uyflux_Fa;
                rhou33.at(y_ind) += m_rho(dof_ind) * uzflux_Fa * uzflux_Fa;
                rhou12.at(y_ind) += m_rho(dof_ind) * uxflux_Fa * uyflux_Fa;
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
        }
    }
    // average of neighboring points if there were no points
    for (auto iy : nonumbers) {
        if (iy == 0) {iy_l = 0; iy_u = iy + 1;}
        else if (iy == m_nofCoordinates - 1) {iy_l = iy - 1; iy_u = iy;}
        else {iy_l = iy - 1; iy_u = iy + 1;}
        dy_u = m_yCoordinates.at(iy_u) - m_yCoordinates.at(iy);
        dy_l = m_yCoordinates.at(iy) - m_yCoordinates.at(iy_l);
        window_size = dy_u + dy_l;
        rhou11.at(iy) = 1 / window_size * (rhou11.at(iy_u) * dy_u + rhou11.at(iy_l) * dy_l);
        rhou22.at(iy) = 1 / window_size * (rhou22.at(iy_u) * dy_u + rhou22.at(iy_l) * dy_l);
        rhou33.at(iy) = 1 / window_size * (rhou33.at(iy_u) * dy_u + rhou33.at(iy_l) * dy_l);
        rhou12.at(iy) = 1 / window_size * (rhou12.at(iy_u) * dy_u + rhou12.at(iy_l) * dy_l);
    }

    // calculate Rij
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        m_R11.at(iy) = rhou11.at(iy) / rho_average.at(iy);
        m_R22.at(iy) = rhou22.at(iy) / rho_average.at(iy);
        m_R33.at(iy) = rhou33.at(iy) / rho_average.at(iy);
        m_R12.at(iy) = rhou12.at(iy) / rho_average.at(iy);
        // calculate turbulent kinetic energy
        m_K.at(iy) = (m_R11.at(iy) + m_R22.at(iy) + m_R33.at(iy)) / 2; // https://www.osti.gov/pages/servlets/purl/1580489
        // calculate anisotropy tensor along y
        b11vec.at(iy) = (m_R11.at(iy) - 2. / 3 * m_K.at(iy) * 1) / (2 * m_K.at(iy));
        b22vec.at(iy) = (m_R22.at(iy) - 2. / 3 * m_K.at(iy) * 1) / (2 * m_K.at(iy));
        b12vec.at(iy) = (m_R12.at(iy) - 2. / 3 * m_K.at(iy) * 0) / (2 * m_K.at(iy));
    }

    if (is_MPI_rank_0()) {
        // calculate y-derivatives
        vector<double> dUdy_abs = derivative(umag_average);
        vector<double> uxFavreDy = derivative(ux_favre);
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            dUdy_abs.at(iy) = abs(dUdy_abs.at(iy));
            growthrate_integrand.at(iy) = rho_average.at(iy) * m_R12.at(iy) * uxFavreDy.at(iy);
        }

        // calculate y-integrals
        double momentumthickness_integral = integrate(momentumthickness_integrand);
        double growthrate_integral = integrate(growthrate_integrand);

        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
            // calculate anisotropy tensor elements
            if (iy == 0) { // left side: trapezoidal rule
                window_size = (m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy));
                m_b11 += window_size * 0.5 * (b11vec.at(iy) + b11vec.at(iy + 1));
                m_b22 += window_size * 0.5 * (b22vec.at(iy) + b22vec.at(iy + 1));
                m_b12 += window_size * 0.5 * (b12vec.at(iy) + b12vec.at(iy + 1));
            } else if (iy == m_nofCoordinates - 1) {
                window_size = (m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1));
                m_b11 += window_size * 0.5 * (b11vec.at(iy) + b11vec.at(iy -1));
                m_b22 += window_size * 0.5 * (b22vec.at(iy) + b22vec.at(iy -1));
                m_b12 += window_size * 0.5 * (b12vec.at(iy) + b12vec.at(iy -1));
            } else { // other: simpson rule
                window_size = 0.5 * abs(m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1));
                m_b11 += window_size * (b11vec.at(iy - 1) + 4 * b11vec.at(iy) + b11vec.at(iy + 1)) / 6;
                m_b22 += window_size * (b22vec.at(iy - 1) + 4 * b22vec.at(iy) + b22vec.at(iy + 1)) / 6;
                m_b12 += window_size * (b12vec.at(iy - 1) + 4 * b12vec.at(iy) + b12vec.at(iy + 1)) / 6;
            }
        }

        // calculate vorticity thickness
        m_currentDeltaOmega = 2 /*dU*/ / *max_element(begin(dUdy_abs), end(dUdy_abs));
        m_ReOmega = 1 /*rho0*/ * 2 /*dU*/ * m_Re0 * m_currentDeltaOmega;

        // calculate momentum thickness
        m_currentDeltaTheta = momentumthickness_integral / (rho0 * pow(dU_calculated, 2));
        m_deltaThetaGrowth = growthrate_integral * (-2) / (rho0 * pow(dU_calculated, 3));

        // output to console
//        cout << "dUdy_abs: ";
//        for (size_t yi = 0; yi < m_nofCoordinates; yi++) {
//            cout << dUdy_abs.at(yi) << ",";
//        } cout << endl;
        cout << "IT: " << m_solver.getIteration()
            << ", t: " << m_solver.getTime()
            << ", dU: " << dU_calculated
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

double ShearLayerStats::integrate(vector<double> integrand) {
    double integral = 0;
    double window_size;
    double fi;
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        if (iy == 0) { // left side: trapezoidal rule
            window_size = (m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy));
            fi = (integrand.at(iy) + integrand.at(iy + 1)) / 2;
        } else if (iy == m_nofCoordinates - 1) {
            window_size = (m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1));
            fi = (integrand.at(iy) + integrand.at(iy - 1)) / 2;
        } else { // other: simpson rule
            window_size = 0.5 * (m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1));
            fi = (integrand.at(iy - 1) + 4 * integrand.at(iy) + integrand.at(iy + 1)) / 6;
        }
        integral += window_size * fi;
    }
    return integral;
}

vector<double> ShearLayerStats::derivative(vector<double> values) {
    vector<double> derivat(m_nofCoordinates, 0);
    double dy;
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        if (iy == 0) { // forward
            dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
            derivat.at(iy) = (values.at(iy + 1) - values.at(iy)) / dy;
        } else if (iy == m_nofCoordinates-1) { // backward
            dy = m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1);
            derivat.at(iy) = (values.at(iy) - values.at(iy - 1)) / dy;
        } else { // other: central
            dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1);
            derivat.at(iy) = (values.at(iy + 1) -  values.at(iy - 1)) / dy;
        }
    }

    return derivat;
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
