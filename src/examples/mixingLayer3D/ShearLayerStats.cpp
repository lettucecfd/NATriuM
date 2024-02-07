/*
 * ShearLayerStats.cpp
 *
 *  Created on: 06.10.2022
 *      Author: dwilde3m
 */

#include "ShearLayerStats.h"

#include "mpi.h"
#include <utility>
#include <filesystem>
#include "deal.II/base/mpi.h"
//#include "natrium/collision_advanced/AuxiliaryCollisionFunctions.h"

namespace natrium {

//ShearLayerStats::ShearLayerStats(CompressibleCFDSolver<3> &solver, std::string outdir, double starting_delta_theta,
//                                 double starting_Re, size_t reflevel, vector<unsigned int> repetitions) :
ShearLayerStats::ShearLayerStats(CFDSolver<3> &solver, std::string outdir, double starting_delta_theta,
                                 double starting_Re, size_t reflevel, vector<unsigned int> repetitions) :
        DataProcessor<3>(solver),
        m_Re0(starting_Re), m_u(solver.getVelocity()), m_rho(solver.getDensity()),
        m_outDir(outdir),
        m_filename(scalaroutfile(solver.getConfiguration()->getOutputDirectory())),
        m_vectorfilename(vectoroutfile(solver.getConfiguration()->getOutputDirectory())),
        m_reflevel(reflevel),
        m_currentDeltaTheta_Fa(starting_delta_theta), m_currentDeltaOmega(0.41),
        m_b11(0), m_b22(0), m_b33(0), m_b12(0), m_b13(0), m_b23(0), m_K_integrated(0), m_reps(repetitions),
        m_no_stats(solver.getConfiguration()->getNoStatsInterval()) {

    TimerOutput::Scope timer_section(Timing::getTimer(), "Shearlayer Reporter");

    m_yCoordsUpToDate = false;
    m_nofCoordinates = 0;

    if (solver.getIterationStart() > 0) {
        if (is_MPI_rank_0()) {
            m_tableFile = boost::make_shared<std::fstream>(m_filename, std::fstream::out | std::fstream::app);
            m_vectorFile = boost::make_shared<std::fstream>(m_vectorfilename, std::fstream::out | std::fstream::app);
        }
    } else {
        if (is_MPI_rank_0()) {
            std::filesystem::path out_dir(solver.getConfiguration()->getOutputDirectory() + "/stats");
            std::filesystem::create_directory(out_dir);
            m_tableFile = boost::make_shared<std::fstream>(m_filename, std::fstream::out);
            m_vectorFile = boost::make_shared<std::fstream>(m_vectorfilename, std::fstream::out);
        }
    }

    roundtol = pow(10, -solver.getConfiguration()->getCoordsRound())
            * CFDSolverUtilities::getMinimumDoFDistanceGLC<3>(*m_solver.getProblemDescription()->getMesh(),
                                                              m_solver.getConfiguration()->getSedgOrderOfFiniteElement());
    calculateDeltas(starting_delta_theta);
    updateYValues();
    testIntegration();
    testDerivate();
    calculateRhoU();
    write_tn();
    write_console();
}

bool ShearLayerStats::isMYCoordsUpToDate() const {
    return m_yCoordsUpToDate;
}

void ShearLayerStats::testIntegration() {
    vector<double> ones(m_nofCoordinates, 1.);
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Integration test of int_-dOmeg^dOmeg(1) returned "
                      << integrate(ones, -m_currentDeltaOmega, m_currentDeltaOmega)
                      << ". Expected " << 2*m_currentDeltaOmega << endl;
        LOG(DETAILED) << "Integration test of int_-dT0^dT0(1) returned "
                      << integrate(ones, -m_currentDeltaTheta_Fa, m_currentDeltaTheta_Fa)
                      << ". Expected " << 2*m_currentDeltaTheta_Fa << endl;
        LOG(DETAILED) << "Integration test of int_-1.222^1.222(1) returned "
                      << integrate(ones, -1.222, 1.222)
                      << ". Expected " << 2*1.222 << endl;
        LOG(DETAILED) << "Integration test of int_-1^1(1) returned "
                      << integrate(ones, -1, 1)
                      << ". Expected 2" << endl;
    }
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Integration test of int(1) returned "
                      << integrate(ones)
                      << ". Expected " << m_ly << endl;
    }
}

void ShearLayerStats::testDerivate() {
    vector<double> f(m_nofCoordinates);
    vector<double> deri(m_nofCoordinates);
    vector<double> diff(m_nofCoordinates);
    size_t iy;
    double fac;
    double ideal_result;
    // constant
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        f.at(iy) = 1;
    }
    deri = derivative(f);
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        ideal_result = 0;
        diff.at(iy) = pow((deri.at(iy) - ideal_result)/ideal_result, 2);
    }
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Derivative of f(x) = 1" << endl
                      << "RMS diff to  f'(x) = 0: "
                      << accumulate(diff.begin(), diff.end(), 0.) << ", expected 0." << endl;
    }
    // linear
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        f.at(iy) = m_yCoordinates.at(iy);
    }
    deri = derivative(f);
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        ideal_result = 1;
        diff.at(iy) = pow((deri.at(iy) - ideal_result)/ideal_result, 2);
    }
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Derivative of f(x) = x" << endl
                      << "RMS diff to  f'(x) = 1: "
                      << accumulate(diff.begin(), diff.end(), 0.) << ", expected 0." << endl;
    }
    fac = 4/m_ly;
    // parabola
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        f.at(iy) = pow(fac * m_yCoordinates.at(iy),2);
    }
    deri = derivative(f);
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        ideal_result = 2*fac*(fac*m_yCoordinates.at(iy));
        diff.at(iy) = pow((deri.at(iy) - ideal_result)/ideal_result, 2);
    }
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Derivative of f(x) = (" << fac << "x)^2" << endl
                      << "RMS diff to  f'(x) = 2*" << fac << "*" << fac << "*x: "
                      << accumulate(diff.begin(), diff.end(), 0.) << ", expected 0." << endl;
    }
    // sine
    fac = 2*3.1415/m_ly;
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        f.at(iy) = sin(fac*m_yCoordinates.at(iy));
    }
    deri = derivative(f);
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        ideal_result = fac*cos(fac*m_yCoordinates.at(iy));
        diff.at(iy) = pow((deri.at(iy) - ideal_result)/ideal_result, 2);
    }
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Derivative of f(x) = sin(" << fac << "*x)" << endl
                      << "RMS diff to  f'(x) = " << fac << "*cos(" << fac << "*x) : "
                      << accumulate(diff.begin(), diff.end(), 0.) << ", expected 0." << endl;
    }
    // cosine
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        f.at(iy) = cos(fac*m_yCoordinates.at(iy));
    }
    deri = derivative(f);
    for (iy = 0; iy < m_nofCoordinates; iy++) {
        ideal_result = fac*(-sin(fac*m_yCoordinates.at(iy)));
        diff.at(iy) = pow((deri.at(iy) - ideal_result)/ideal_result, 2);
    }
    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Derivative of f(x) = cos(" << fac << "*x)" << endl
                      << "RMS diff to  f'(x) = " << fac << "*-sin(" << fac << "*x) : "
                      << accumulate(diff.begin(), diff.end(), 0.) << ", expected 0." << endl;
    }
    bool debug = false;
    if (debug) {
        LOG(DETAILED) << "y:" << endl << "[";
        for (iy = 0; iy < m_nofCoordinates-1; iy++) {
            LOG(DETAILED) << m_yCoordinates.at(iy) << ", ";
        }
        LOG(DETAILED) << m_yCoordinates.at(m_nofCoordinates-1) << "]" << endl;
        LOG(DETAILED) << "cos(y):" << endl << "[";
        for (iy = 0; iy < m_nofCoordinates-1; iy++) {
            LOG(DETAILED) << f.at(iy) << ", ";
        }
        LOG(DETAILED) << f.at(m_nofCoordinates-1) << "]" << endl;
        LOG(DETAILED) << "(cos(y))':" << endl << "[";
        for (iy = 0; iy < m_nofCoordinates-1; iy++) {
            LOG(DETAILED) << deri.at(iy) << ", ";
        }
        LOG(DETAILED) << deri.at(m_nofCoordinates-1) << "]" << endl;
    }
}

void ShearLayerStats::calculateDeltas(double dT0) {
//    size_t itest = 0;
    for (size_t dim = 0; dim < 3; ++dim) {
        boost::shared_ptr<AdvectionOperator<3> > advection = m_solver.getAdvectionOperator();

        //////////////////////////
        // Calculate coordinates ////
        //////////////////////////
        std::set<double, own_double_less> coords;

        // get y coordinates
        const dealii::UpdateFlags update_flags = dealii::update_quadrature_points;
        const dealii::DoFHandler<3> &dof_handler = *(advection->getDoFHandler());
        dealii::FEValues<3> fe_values(advection->getMapping(), *(advection->getFe()),
                                      advection->getSupportPointEvaluation(), update_flags);
        typename dealii::DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                fe_values.reinit(cell);
                const std::vector<dealii::Point<3> > &quad_points = fe_values.get_quadrature_points();
                for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
//                    coords.insert(floor(quad_points.at(i)(dim) * (nround)) / (nround));
                    coords.insert(quad_points.at(i)(dim));
                }
            }
        }
        //communicate list of y values
        size_t n_coords = coords.size();
        size_t max_ncoords = dealii::Utilities::MPI::min_max_avg(n_coords, MPI_COMM_WORLD).max;
        auto *sendbuf = (double *) malloc(max_ncoords * sizeof(double));

        // fill send buffer with y coordinates
        size_t i = 0;
        for (double coord: coords) {
            sendbuf[i] = coord;
            i++;
        }
        for (; i < max_ncoords; i++) {
            sendbuf[i] = *coords.begin();
        }
        // allocate read buffer
        double *recvbuf = (double *) malloc(
                dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) * max_ncoords * sizeof(double));
        // allgather: distribute all ycoords to all processors
        MPI_Allgather(sendbuf, max_ncoords, MPI_DOUBLE, recvbuf, max_ncoords, MPI_DOUBLE, MPI_COMM_WORLD);

        // remove double entries
        std::set<double> coords_gathered = {};
        for (i = 0; i < max_ncoords * dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); i++) {
            bool is_unique = true;
            for (auto coord : coords_gathered) {
                if (std::abs(recvbuf[i] - coord) < roundtol) {
                    is_unique = false;
                    break;
                }
            } // all so-far-unique coordinates
            if (is_unique) {
                m_nofCoordinates_all.at(dim) += 1;
//                coords_gathered.resize(m_nofCoordinates_all.at(dim));
                coords_gathered.insert(recvbuf[i]);
            }
        } // all gathered coordinates

        // fill member variables
        i = 0;
        for (double it: coords_gathered) {
            if (dim == 0) {
                m_xCoordinates.push_back(it);
                m_xCoordinateToIndex.insert(std::make_pair(it, i));
            }
            if (dim == 1) {
                m_yCoordinates.push_back(it);
                m_yCoordinateToIndex.insert(std::make_pair(it, i));
            }
            if (dim == 2) {
                m_zCoordinates.push_back(it);
                m_zCoordinateToIndex.insert(std::make_pair(it, i));
            }
            i++;
        }
        // free
        free(sendbuf);
        free(recvbuf);
        // finished
    }
    vector<double> mindeltas(3,10000), maxdeltas(3,0); // double mindx=0, maxdx=0, mindy=0, maxdy=0, mindz=0, maxdz=0;
    vector<vector<double>*> all_coords;
    all_coords.push_back(&m_xCoordinates); all_coords.push_back(&m_yCoordinates); all_coords.push_back(&m_zCoordinates);
    vector<double> deltas;

    for (size_t dim = 0; dim < 3; ++dim) {
        //// calculate deltas from coordinates
        deltas.resize(m_nofCoordinates_all.at(dim) - 1);
        for (size_t i = 0; i < m_nofCoordinates_all.at(dim)-1; ++i) {
            deltas.at(i) = all_coords.at(dim)->at(i+1) - all_coords.at(dim)->at(i);
        }
        mindeltas.at(dim) = *std::min_element(deltas.begin(), deltas.end());
        maxdeltas.at(dim) = *std::max_element(deltas.begin(), deltas.end());
        // communicate
        mindeltas.at(dim) = dealii::Utilities::MPI::min_max_avg(mindeltas.at(dim), MPI_COMM_WORLD).min;
        maxdeltas.at(dim) = dealii::Utilities::MPI::min_max_avg(maxdeltas.at(dim), MPI_COMM_WORLD).max;
    }
    double dofmin, dofmax;
    if (m_solver.getConfiguration()->getSupportPoints() == GAUSS_LOBATTO_CHEBYSHEV_POINTS) {
        dofmin = CFDSolverUtilities::getMinimumDoFDistanceGLC<3>(*m_solver.getProblemDescription()->getMesh(),
                                                                        m_solver.getConfiguration()->getSedgOrderOfFiniteElement());
        dofmax = CFDSolverUtilities::getMaximumDoFDistanceGLC<3>(*m_solver.getProblemDescription()->getMesh(),
                                                                        m_solver.getConfiguration()->getSedgOrderOfFiniteElement());
    } else if (m_solver.getConfiguration()->getSupportPoints() == EQUIDISTANT_POINTS) {
        dofmin = CFDSolverUtilities::getMinimumVertexDistance<3>(*m_solver.getProblemDescription()->getMesh());
        dofmax = CFDSolverUtilities::getMinimumVertexDistance<3>(*m_solver.getProblemDescription()->getMesh());
    } else {
        dofmin = 0; dofmax = 0;
    }
    vector<double> mindeltas_verteces = CFDSolverUtilities::getMinimumVertexDistanceDirs<3>(*m_solver.getProblemDescription()->getMesh());
    vector<double> maxdeltas_verteces = CFDSolverUtilities::getMaximumVertexDistanceDirs<3>(*m_solver.getProblemDescription()->getMesh());

    unsigned int nx = m_reps.at(0)*pow(2,m_reflevel), ny = m_reps.at(1)*pow(2,m_reflevel), nz = m_reps.at(2)*pow(2,m_reflevel);
    unsigned int nxp = m_xCoordinates.size(), nyp = m_yCoordinates.size(), nzp = m_zCoordinates.size();
    if (is_MPI_rank_0()) {
        double xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = *std::min_element(m_xCoordinates.begin(), m_xCoordinates.end());
        ymin = *std::min_element(m_yCoordinates.begin(), m_yCoordinates.end());
        zmin = *std::min_element(m_zCoordinates.begin(), m_zCoordinates.end());
        xmax = *std::max_element(m_xCoordinates.begin(), m_xCoordinates.end());
        ymax = *std::max_element(m_yCoordinates.begin(), m_yCoordinates.end());
        zmax = *std::max_element(m_zCoordinates.begin(), m_zCoordinates.end());
        m_lx = xmax - xmin; m_ly = ymax - ymin; m_lz = zmax - zmin;
        LOG(DETAILED) << "::::::---------------------------------------" << endl
                      << "Mesh info after transform() (rounded coordinates to " << roundtol << " for uniqueness): " << endl
                          << " dx in [" << mindeltas_verteces.at(0) << "," << maxdeltas_verteces.at(0) << "], " << endl
                          << " dy in [" << mindeltas_verteces.at(1) << "," << maxdeltas_verteces.at(1) << "], " << endl
                          << " dz in [" << mindeltas_verteces.at(2) << "," << maxdeltas_verteces.at(2) << "]." << endl
                          << " x in [" << xmin << "," << xmax << "], " << endl
                          << " y in [" << ymin << "," << ymax << "], " << endl
                          << " z in [" << zmin << "," << zmax << "]." << endl
                          << " ncells = " << nx << "x" << ny << "x" << nz << " = " << nx*ny*nz << endl
                      << "Integration point distances: " << endl
                          << " dx in [" << mindeltas.at(0) << "," << maxdeltas.at(0) << "], " << endl
                          << " dy in [" << mindeltas.at(1) << "," << maxdeltas.at(1) << "], " << endl
                          << " dz in [" << mindeltas.at(2) << "," << maxdeltas.at(2) << "]." << endl
                          << " DOF distance in [" << dofmin << "," << dofmax << "]." << endl
                          << " npoints = " << nxp << "x" << nyp << "x" << nzp << " = " << nxp*nyp*nzp << endl
                      << "normalized by deltaTheta0: " << endl
                          << " dx in [" << mindeltas.at(0) / dT0 << "," << maxdeltas.at(0) / dT0 << "], " << endl
                          << " dy in [" << mindeltas.at(1) / dT0 << "," << maxdeltas.at(1) / dT0 << "], " << endl
                          << " dz in [" << mindeltas.at(2) / dT0 << "," << maxdeltas.at(2) / dT0 << "]." << endl
                          << " DOF distance in [" << dofmin / dT0 << "," << dofmax / dT0 << "]." << endl
                          << " x in [" << xmin / dT0 << "," << xmax / dT0 << "], " << endl
                          << " y in [" << ymin / dT0 << "," << ymax / dT0 << "], " << endl
                          << " z in [" << zmin / dT0 << "," << zmax / dT0 << "]." << endl
                          << " lx = " << m_lx / dT0 << ", ly = " << m_ly / dT0 << ", lz = " << m_lz / dT0 << endl
                      << "---------------------------------------" << endl;
    }
    bool debug = false;
    if (is_MPI_rank_0() and debug) {
        for (int scaleto = 2; scaleto < 7; scaleto ++) {
            double fac = pow(2., int(m_reflevel) - scaleto);
            LOG(DETAILED) << "::::::---------------------------------------" << endl
                          << "Mesh info for ref-level " << scaleto << " (this was " << m_reflevel
                          << ", so multiplying by " << fac << "): " << endl
                          << " dx in [" << mindeltas_verteces.at(0) * fac << "," << maxdeltas_verteces.at(0) * fac
                          << "], " << endl
                          << " dy in [" << mindeltas_verteces.at(1) * fac << "," << maxdeltas_verteces.at(1) * fac
                          << "], " << endl
                          << " dz in [" << mindeltas_verteces.at(2) * fac << "," << maxdeltas_verteces.at(2) * fac
                          << "]." << endl
                          << " ncells = " << int(nx/fac) << "x" << int(ny/fac) << "x" << int(nz/fac) << " = " << int(nx*ny*nz/pow(fac,3)) << endl
                          << "Integration point distances: " << endl
                          << " dx in [" << mindeltas.at(0) * fac << "," << maxdeltas.at(0) * fac << "], " << endl
                          << " dy in [" << mindeltas.at(1) * fac << "," << maxdeltas.at(1) * fac << "], " << endl
                          << " dz in [" << mindeltas.at(2) * fac << "," << maxdeltas.at(2) * fac << "]." << endl
                          << " DOF distance in [" << dofmin * fac << "," << dofmax * fac << "]." << endl
                          << " npoints = " << int(nxp/fac) << "x" << int(nyp/fac) << "x" << int(nzp/fac) << " = " << int(nxp*nyp*nzp/pow(fac,3)) << endl
                          << "normalized by deltaTheta0: " << endl
                          << " dx in [" << mindeltas.at(0) * fac / dT0 << "," << maxdeltas.at(0) * fac / dT0 << "], "
                          << endl
                          << " dy in [" << mindeltas.at(1) * fac / dT0 << "," << maxdeltas.at(1) * fac / dT0 << "], "
                          << endl
                          << " dz in [" << mindeltas.at(2) * fac / dT0 << "," << maxdeltas.at(2) * fac / dT0 << "]."
                          << endl
                          << " DOF distance in [" << dofmin * fac / dT0 << "," << dofmax * fac / dT0 << "]." << endl
                          << "---------------------------------------" << endl;
        }
    }
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
//                y_coords.insert(floor(quad_points.at(i)(1)*nround)/nround);
                y_coords.insert(quad_points.at(i)(1));
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
    std::set<double> y_coords_gathered = {};
    for (i = 0; i < max_nycoords * dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); i++) {
        bool is_unique = true;
        for (auto coord : y_coords_gathered) {
            if (std::abs(recvbuf[i] - coord) < roundtol) {
                is_unique = false;
                break;
            }
        } // all so-far-unique coordinates
        if (is_unique) {
            m_nofCoordinates += 1;
//                coords_gathered.resize(m_nofCoordinates_all.at(dim));
            y_coords_gathered.insert(recvbuf[i]);
        }
    } // all gathered coordinates

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
    TimerOutput::Scope timer_section(Timing::getTimer(), "Shearlayer Reporter");
    size_t iteration = m_solver.getIteration();
    if (int(iteration) < m_no_stats) {
        if (is_MPI_rank_0()) {
            LOG(DETAILED) << "No stats at " << iteration << ", only after " << m_no_stats << endl;
        }
    }
    if (((iteration == 1) or (iteration == 10) or (iteration == 50) or (iteration == 100) or (iteration == 500)
        or (iteration % 1000 == 0)) or ((iteration % m_solver.getConfiguration()->getOutputShearLayerInterval() == 0)
        and (int(iteration) > m_no_stats)))
        {
            if (!isMYCoordsUpToDate()) {
                updateYValues();
            }
            calculateRhoU();
            write_tn();
            write_console();
        }
}

void ShearLayerStats::calculateRhoU() {
    //////////////////////////
    // Calculate averages ////
    //////////////////////////
    // initialize vectors for averages along x and z
    vector<double> rhoux_Re(m_nofCoordinates, 0), rhouy_Re(m_nofCoordinates, 0), rhouz_Re(m_nofCoordinates, 0);
    // initialize vectors for combined sizes
    vector<double> rhoU11Flux(m_nofCoordinates, 0), rhoU22Flux(m_nofCoordinates, 0), rhoU33Flux(m_nofCoordinates, 0),
    rhoU12Flux(m_nofCoordinates, 0), rhoU13Flux(m_nofCoordinates, 0), rhoU23Flux(m_nofCoordinates, 0),
    b11vec(m_nofCoordinates,0), b22vec(m_nofCoordinates,0), b33vec(m_nofCoordinates,0),
    b12vec(m_nofCoordinates,0), b13vec(m_nofCoordinates,0), b23vec(m_nofCoordinates,0),
    growthrate_integrand(m_nofCoordinates,0);
    vector<vector<double>*> rhoUijFlux_set, bijvec_set, Rij_set;
    rhoUijFlux_set.push_back(&rhoU11Flux); rhoUijFlux_set.push_back(&rhoU22Flux); rhoUijFlux_set.push_back(&rhoU33Flux);
    rhoUijFlux_set.push_back(&rhoU12Flux); rhoUijFlux_set.push_back(&rhoU13Flux); rhoUijFlux_set.push_back(&rhoU23Flux);
    bijvec_set.push_back(&b11vec); bijvec_set.push_back(&b22vec); bijvec_set.push_back(&b33vec);
    bijvec_set.push_back(&b12vec); bijvec_set.push_back(&b13vec); bijvec_set.push_back(&b23vec);
    Rij_set.push_back(&m_R11); Rij_set.push_back(&m_R22); Rij_set.push_back(&m_R33);
    Rij_set.push_back(&m_R12); Rij_set.push_back(&m_R13); Rij_set.push_back(&m_R23);
    vector<double*> bij_set;
    bij_set.push_back(&m_b11); bij_set.push_back(&m_b22); bij_set.push_back(&m_b33);
    bij_set.push_back(&m_b12); bij_set.push_back(&m_b13); bij_set.push_back(&m_b23);
    vector<vector<double>*> u_Fa_set, rhou_Re_set;
    u_Fa_set.push_back(&ux_Fa); u_Fa_set.push_back(&uy_Fa); u_Fa_set.push_back(&uz_Fa);
    rhou_Re_set.push_back(&rhoux_Re); rhou_Re_set.push_back(&rhouy_Re); rhou_Re_set.push_back(&rhouz_Re);

    // resize to fit length of y and set elements to 0
    for (size_t dim = 0; dim < 3; dim++) {
        u_Fa_set.at(dim)->resize(m_nofCoordinates); std::fill(u_Fa_set.at(dim)->begin(), u_Fa_set.at(dim)->end(), 0);
        rhou_Re_set.at(dim)->resize(m_nofCoordinates); std::fill(rhou_Re_set.at(dim)->begin(), rhou_Re_set.at(dim)->end(), 0);
    }
    m_number.resize(m_nofCoordinates); std::fill(m_number.begin(), m_number.end(), 0);
    umag_Re.resize(m_nofCoordinates); std::fill(umag_Re.begin(), umag_Re.end(), 0);
    ux_Re.resize(m_nofCoordinates); std::fill(ux_Re.begin(), ux_Re.end(), 0);
    rho_Re.resize(m_nofCoordinates); std::fill(rho_Re.begin(), rho_Re.end(), 0);
    for (size_t ij = 0; ij < 6; ij++) {
        Rij_set.at(ij)->resize(m_nofCoordinates);
        std::fill(Rij_set.at(ij)->begin(), Rij_set.at(ij)->end(), 0);
    }

    m_K.resize(m_nofCoordinates);
    momentumthickness_integrand.resize(m_nofCoordinates);
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
//                y = floor(quad_points.at(i)(1)*nround)/nround;
                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                m_number.at(y_ind) += 1;
                // fill value vector
                for (size_t dim = 0; dim < 3; dim++) {
                    rhou_Re_set.at(dim)->at(y_ind) += m_rho(dof_ind) * m_u.at(dim)(dof_ind);
                }
                rho_Re.at(y_ind) += m_rho(dof_ind);
                ux_Re.at(y_ind) += m_u.at(0)(dof_ind);
                umag_Re.at(y_ind) += sqrt(pow(m_u.at(0)(dof_ind), 2) + pow(m_u.at(1)(dof_ind), 2) + pow(m_u.at(2)(dof_ind), 2));
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */

    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // add n_points(iy), rho*ux(iy), rho(iy), etc. from different MPI processes
        m_number.at(iy) = dealii::Utilities::MPI::sum(m_number.at(iy), MPI_COMM_WORLD);
        // average over number of points at y
        if (m_number.at(iy) != 0) {
            for (size_t dim = 0; dim < 3; dim++) {
                rhou_Re_set.at(dim)->at(iy) = dealii::Utilities::MPI::sum(rhou_Re_set.at(dim)->at(iy), MPI_COMM_WORLD);
                rhou_Re_set.at(dim)->at(iy) /= m_number.at(iy);
            }
            rho_Re.at(iy) = dealii::Utilities::MPI::sum(rho_Re.at(iy), MPI_COMM_WORLD);
            umag_Re.at(iy) = dealii::Utilities::MPI::sum(umag_Re.at(iy), MPI_COMM_WORLD);
            ux_Re.at(iy) = dealii::Utilities::MPI::sum(ux_Re.at(iy), MPI_COMM_WORLD);
            rho_Re.at(iy) /= m_number.at(iy);
            umag_Re.at(iy) /= m_number.at(iy);
            ux_Re.at(iy) /= m_number.at(iy);
        } else nonumbers.push_back(iy);
    }
    // average of neighboring points if there were no points
    size_t iy_u, iy_l;
    float dy_u, dy_l;
    float window_size;
    for (auto iy : nonumbers) {
        if (iy == 0) {iy_l = 0; iy_u = iy + 1;}
        else if (iy == m_nofCoordinates - 1) {iy_l = iy - 1; iy_u = iy;}
        else {iy_l = iy - 1; iy_u = iy + 1;}
        dy_u = m_yCoordinates.at(iy_u) - m_yCoordinates.at(iy);
        dy_l = m_yCoordinates.at(iy) - m_yCoordinates.at(iy_l);
        window_size = dy_u + dy_l;
        for (size_t dim = 0; dim < 3; dim++) {
            rhou_Re_set.at(dim)->at(iy) = 1 / window_size * (rhou_Re_set.at(dim)->at(iy_u) * dy_u + rhou_Re_set.at(dim)->at(iy_l) * dy_l);
        }
        rho_Re.at(iy) = 1 / window_size * (rho_Re.at(iy_u) * dy_u + rho_Re.at(iy_l) * dy_l);
        umag_Re.at(iy) = 1 / window_size * (umag_Re.at(iy_u) * dy_u + umag_Re.at(iy_l) * dy_l);
        ux_Re.at(iy) = 1 / window_size * (ux_Re.at(iy_u) * dy_u + ux_Re.at(iy_l) * dy_l);
    }

    // ui_favre and momentumthickness_integrand
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        for (size_t dim = 0; dim < 3; dim++) {
            u_Fa_set.at(dim)->at(iy) = rhou_Re_set.at(dim)->at(iy) / rho_Re.at(iy);
        }
    }
    auto [minUx, maxUx] = std::minmax_element(begin(ux_Fa), end(ux_Fa));
    m_dUx = *maxUx - *minUx;
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        momentumthickness_integrand.at(iy) = rho_Re.at(iy) * (m_dU0 / 2 - ux_Fa.at(iy)) * (m_dU0 / 2 + ux_Fa.at(iy));
    }

    // calculate flux of ux, uy, uz at all points
    vector<double> uflux_Fa_set(4, 0);
    for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<3> > &quad_points = fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
//                y = floor(quad_points.at(i)(1)*nround)/nround;
                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);
                dof_ind = local_indices.at(i);
                // calculate fluxes
                for (size_t dim = 0; dim < 3; dim++) {
                    uflux_Fa_set.at(dim) = m_u.at(dim)(dof_ind) - u_Fa_set.at(dim)->at(y_ind);
                }
                // fill value vector
                for (size_t ij = 0; ij < 3; ij++) {
                    rhoUijFlux_set.at(ij)->at(y_ind) += m_rho(dof_ind) * uflux_Fa_set.at(ij) * uflux_Fa_set.at(ij);
                }
                rhoUijFlux_set.at(3)->at(y_ind) += m_rho(dof_ind) * uflux_Fa_set.at(0) * uflux_Fa_set.at(1); // rhou1'u2'
                rhoUijFlux_set.at(4)->at(y_ind) += m_rho(dof_ind) * uflux_Fa_set.at(0) * uflux_Fa_set.at(2); // rhou1'u3'
                rhoUijFlux_set.at(5)->at(y_ind) += m_rho(dof_ind) * uflux_Fa_set.at(1) * uflux_Fa_set.at(2); // rhou2'u3'
            } /* for all quadrature points */
        } /* if locally owned */
    } /* for all cells */
    // communicate
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        // average over number of points at y
        if (m_number.at(iy) != 0) {
            for (size_t ij = 0; ij < 6; ij++) {
                rhoUijFlux_set.at(ij)->at(iy) = dealii::Utilities::MPI::sum(rhoUijFlux_set.at(ij)->at(iy), MPI_COMM_WORLD);
                rhoUijFlux_set.at(ij)->at(iy) /= m_number.at(iy);
            }
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
        for (size_t ij = 0; ij < 6; ij++) {
            rhoUijFlux_set.at(ij)->at(iy) = 1 / window_size * (rhoUijFlux_set.at(ij)->at(iy_u) * dy_u + rhoUijFlux_set.at(ij)->at(iy_l) * dy_l);
        }
    }

    // calculate Rij
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        for (size_t ij = 0; ij < 6; ij++) {
            Rij_set.at(ij)->at(iy) = rhoUijFlux_set.at(ij)->at(iy) / rho_Re.at(iy);
        }
        // calculate turbulent kinetic energy
        m_K.at(iy) = (m_R11.at(iy) + m_R22.at(iy) + m_R33.at(iy)) / 2; // https://www.osti.gov/pages/servlets/purl/1580489
        // calculate anisotropy tensor along y
        for (size_t ij = 0; ij < 3; ij++) {
            bijvec_set.at(ij)->at(iy) = (Rij_set.at(ij)->at(iy) - 2. / 3 * m_K.at(iy)) / (2 * m_K.at(iy));
        }
        bijvec_set.at(3)->at(iy) = Rij_set.at(3)->at(iy) / (2 * m_K.at(iy));
        bijvec_set.at(4)->at(iy) = Rij_set.at(4)->at(iy) / (2 * m_K.at(iy));
        bijvec_set.at(5)->at(iy) = Rij_set.at(5)->at(iy) / (2 * m_K.at(iy));
    }

    if (is_MPI_rank_0()) {
        // calculate y-derivatives
//        vector<double> dUdy_abs = derivative(umag_Re);
        vector<double> ux_Fa_dy = derivative(ux_Fa);
        for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
//            dUdy_abs.at(iy) = abs(dUdy_abs.at(iy));
            growthrate_integrand.at(iy) = rho_Re.at(iy) * m_R12.at(iy) * ux_Fa_dy.at(iy);
            ux_Fa_dy.at(iy) = abs(ux_Fa_dy.at(iy));
        }

        // calculate y-integrals
        double momentumthickness_integral_Fa = integrate(momentumthickness_integrand);
        double growthrate_integral = integrate(growthrate_integrand);

        // calculate vorticity thickness
        m_currentDeltaOmega = m_dU0 /*dU*/ / *max_element(begin(ux_Fa_dy), end(ux_Fa_dy));
        m_ReOmega = rho0 * m_Re0 * m_dU0 * m_currentDeltaOmega;

        // calculate momentum thickness
        m_currentDeltaTheta_Fa = momentumthickness_integral_Fa / (rho0 * pow(m_dU0, 2));
        m_deltaThetaGrowth = growthrate_integral * (-2) / (rho0 * pow(m_dU0, 3));

        // integrate bij over vorticity thickness
        for (size_t ij = 0; ij < 6; ij++) {
            *bij_set.at(ij) = 1 / (2*m_currentDeltaOmega)
                    * integrate(*bijvec_set.at(ij), -m_currentDeltaOmega, m_currentDeltaOmega);
        }
        m_K_integrated = 1 / (2*m_currentDeltaOmega)
                         * integrate(m_K, -m_currentDeltaOmega, m_currentDeltaOmega);
    }
    auto [minR11_pos, maxR11_pos] = std::minmax_element(begin(m_R11), end(m_R11));
    min_R11 = *minR11_pos; max_R11 = *maxR11_pos;
    auto [minR22_pos, maxR22_pos] = std::minmax_element(begin(m_R22), end(m_R22));
    min_R22 = *minR22_pos; max_R22 = *maxR22_pos;
    auto [minR33_pos, maxR33_pos] = std::minmax_element(begin(m_R33), end(m_R33));
    min_R33 = *minR33_pos; max_R33 = *maxR33_pos;
    auto [minR12_pos, maxR12_pos] = std::minmax_element(begin(m_R12), end(m_R12));
    min_R12 = *minR12_pos; max_R12 = *maxR12_pos;
    auto [minR13_pos, maxR13_pos] = std::minmax_element(begin(m_R13), end(m_R13));
    min_R13 = *minR13_pos; max_R13 = *maxR13_pos;
    auto [minR23_pos, maxR23_pos] = std::minmax_element(begin(m_R23), end(m_R23));
    min_R23 = *minR23_pos; max_R23 = *maxR23_pos;
}

double ShearLayerStats::integrate(vector<double> integrand) {
    double integral = 0;
    double window_size;
    double fi;
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        if (iy == 0) { // left side: trapezoidal rule
            window_size = (m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy)) / 2;
            fi = (integrand.at(iy) + integrand.at(iy + 1)) / 2;
        } else if (iy == m_nofCoordinates - 1) {
            window_size = (m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1)) / 2;
            fi = (integrand.at(iy) + integrand.at(iy - 1)) / 2;
        } else { // other: simpson rule
            window_size = 0.5 * (m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy - 1));
            fi = (integrand.at(iy - 1) + 4 * integrand.at(iy) + integrand.at(iy + 1)) / 6;
        }
        integral += window_size * fi;
    }
    return integral;
}

double ShearLayerStats::integrate(vector<double> integrand, double ymin, double ymax) {
    double integral = 0;
    double window_size, fi, yi, y_upper2, y_lower2, y_upper, y_lower;
    //// TODO: additional extrapolation
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        yi = m_yCoordinates.at(iy);
        window_size = 0;
        fi = 0;
        if (iy == 0) {
            if ((ymin < yi) and (yi < ymax)) {// left side: trapezoidal rule
                window_size = m_yCoordinates.at(iy + 1) - yi;
                fi = (integrand.at(iy) + integrand.at(iy + 1)) / 2;
            }
        } else if (iy == m_nofCoordinates - 1) {
            if ((ymin < yi) and (yi < ymax)) {
                window_size = yi - m_yCoordinates.at(iy - 1);
                fi = (integrand.at(iy) + integrand.at(iy - 1)) / 2;
            }
        } else if ((iy > 0) and (iy < m_nofCoordinates - 1)){
            y_upper = m_yCoordinates.at(iy + 1);
            y_lower = m_yCoordinates.at(iy - 1);
            if ((ymin < y_upper) and (y_lower < ymax)) { // other: simpson rule
                y_upper2 = (y_upper + yi) / 2;
                y_lower2 = (y_lower + yi) / 2;
                fi = (integrand.at(iy - 1) + 4 * integrand.at(iy) + integrand.at(iy + 1)) / 6;
                window_size = std::max(std::min(y_upper2, ymax) - std::max(y_lower2, ymin), 0.);
            }
        }
//        if (is_MPI_rank_0()) {
//            cout << "yi="<<yi<<",ymin="<<ymin<<",ymax="<<ymax<<",yl2="<<y_lower2<<",yu2="<<y_upper2<<",l="<<window_size<<endl;
//        }
        integral += window_size * fi;
    }
    return integral;
}

vector<double> ShearLayerStats::derivative(vector<double> values) {
    vector<double> derivat(m_nofCoordinates, 0);
    double dy;
    double hl;
    double hu;
    for (size_t iy = 0; iy < m_nofCoordinates; iy++) {
        if (iy == 0) { // forward
            dy = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
            derivat.at(iy) = (values.at(iy + 1) - values.at(iy)) / dy;
        } else if (iy == m_nofCoordinates-1) { // backward
            dy = m_yCoordinates.at(iy) - m_yCoordinates.at(iy-1);
            derivat.at(iy) = (values.at(iy) - values.at(iy - 1)) / dy;
        } else { // other: central
//            double h = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
            hl = m_yCoordinates.at(iy) - m_yCoordinates.at(iy - 1);
            hu = m_yCoordinates.at(iy + 1) - m_yCoordinates.at(iy);
            derivat.at(iy) = (values.at(iy + 1)*hl*hl + (hu*hu - hl*hl)*values.at(iy) - hu*hu*values.at(iy - 1)) / (hl*hu*(hl+hu));  // numpy documentation
        }
    }
    return derivat;
}

void ShearLayerStats::write_console() {
    if (is_MPI_rank_0()) {
        std::stringstream log;
//        log.precision(6);
        log << "IT: " << m_solver.getIteration()
            << ", t: " << m_solver.getTime()
            << ", delta_Theta: " << m_currentDeltaTheta_Fa
            << ", growth_rate: " << m_deltaThetaGrowth
            << ", delta_Omega: " << m_currentDeltaOmega
            << ", Re_Omega: " << m_ReOmega
            << ", b11: " << m_b11
            << ", b22: " << m_b22
//            << ", b33: " << m_b33
            << ", b12: " << m_b12;
//            << ", b13: " << m_b13
//            << ", b23: " << m_b23
        LOG(DETAILED) << log.str() << endl;
    }
}

void ShearLayerStats::write_tn() {
    if (is_MPI_rank_0()) {
        string dir = m_solver.getConfiguration()->getOutputDirectory();
        boost::filesystem::path out_dir(dir);
        string filename = "stats/shearlayer_t" + std::to_string(m_solver.getIteration()) + ".txt";
        boost::filesystem::path out_file = out_dir / filename;
        std::ofstream ofs;
        ofs.open(out_file, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
        const string& tnFilename = out_file.string();
        boost::shared_ptr<std::fstream> tnFile = boost::make_shared<std::fstream>(tnFilename, std::fstream::out | std::fstream::app);
        *tnFile << "it t deltaTheta_Fa deltaThetaDot deltaOmega b11 b22 b33 b12 b13 b23 "
                   "R11min R11max R22min R22max R33min R33max "
                   "R12min R12max R13min R13max R23min R23max k_integrated dUx ReOmega" << endl;
        *tnFile << this->m_solver.getIteration() << " " << m_solver.getTime() << " "
                << m_currentDeltaTheta_Fa << " " << m_deltaThetaGrowth << " " << m_currentDeltaOmega
                << " " << m_b11 << " " << m_b22 << " " << m_b33 << " " << m_b12 << " " << m_b13 << " " << m_b23
                << " " << min_R11 << " " << max_R11 << " " << min_R22 << " " << max_R22
                << " " << min_R33 << " " << max_R33 << " " << min_R12 << " " << max_R12
                << " " << min_R13 << " " << max_R13 << " " << min_R23 << " " << max_R23 << " " << m_K_integrated
                << " " << m_dUx << " " << m_ReOmega
                << endl;
        *tnFile << "y: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_yCoordinates.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "ux_Fa: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << ux_Fa.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "rho_Re: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << rho_Re.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "umag_Re: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << umag_Re.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "R11: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_R11.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "R22: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_R22.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "R33: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_R33.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "R12: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_R12.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "R13: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_R13.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "R23: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_R23.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "uy_Fa: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << uy_Fa.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "uz_Fa: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << uz_Fa.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "K: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << m_K.at(iy) << " ";
        } *tnFile << endl;
        *tnFile << "ux_Re: ";
        for (size_t iy = 0; iy < m_nofCoordinates-1; iy++) {
            *tnFile << ux_Re.at(iy) << " ";
        } *tnFile << endl;
    }
}

ShearLayerStats::~ShearLayerStats() = default;

} /* namespace natrium */
