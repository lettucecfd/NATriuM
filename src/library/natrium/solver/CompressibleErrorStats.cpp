//
// Created by dwilde3m on 25.02.20.
//

#include "CompressibleErrorStats.h"
/*
 * CompressibleErrorStats.cpp
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#include "CompressibleErrorStats.h"

#include "BenchmarkCompressibleCFDSolver.h"
#include "deal.II/numerics/vector_tools.h"
#include "deal.II/base/mpi.h"
#include "math.h"

namespace natrium {

    template<size_t dim>
    CompressibleErrorStats<dim>::CompressibleErrorStats(BenchmarkCompressibleCFDSolver<dim> * cfdsolver,
                                const std::string tableFileName) :
            m_solver(cfdsolver), m_filename(tableFileName), m_outputOff(
            tableFileName == "") {
        // set information
        m_iterationNumber = 100000000037;
        m_time = 0.0;
        m_maxVelocityError = 0.0;
        m_maxDensityError = 0.0;
        m_l2VelocityError = 0.0;
        m_l2DensityError = 0.0;
        m_maxUAnalytic = 0.0;

        // create file (if necessary)
        if (m_solver->getIterationStart() > 0) {
            m_errorsTableFile = boost::make_shared<std::fstream>(tableFileName,
                                                                 std::fstream::out | std::fstream::app);
        } else {
            m_errorsTableFile = boost::make_shared<std::fstream>(tableFileName,
                                                                 std::fstream::out);
            printHeaderLine();
        }
        MPI_sync();
    }
    template CompressibleErrorStats<2>::CompressibleErrorStats(BenchmarkCompressibleCFDSolver<2> * cfdsolver,
                                       const std::string tableFileName);
    template CompressibleErrorStats<3>::CompressibleErrorStats(BenchmarkCompressibleCFDSolver<3> * cfdsolver,
                                       const std::string tableFileName);

    template<size_t dim>
    void CompressibleErrorStats<dim>::update() {
        // this function must not be called more often than once per iteration
        // as the data for the analytic solution is constantly overwritten
        // therefor check a marker value that is set by this function (see below)
        if (m_iterationNumber == m_solver->getIteration()) {
            return;
        }


        vector<distributed_vector> num_u;
        distributed_vector num_rho;
        vector<distributed_vector> ana_u;
        distributed_vector ana_rho;
        CFDSolverUtilities::getWriteableDensity(num_rho, m_solver->getDensity(),
                                                m_solver->getAdvectionOperator()->getLocallyOwnedDofs());
        CFDSolverUtilities::getWriteableVelocity(num_u, m_solver->getVelocity(),
                                                 m_solver->getAdvectionOperator()->getLocallyOwnedDofs());
        CFDSolverUtilities::getWriteableDensity(ana_rho, m_solver->getDensity(),
                                                m_solver->getAdvectionOperator()->getLocallyOwnedDofs());
        CFDSolverUtilities::getWriteableVelocity(ana_u, m_solver->getVelocity(),
                                                 m_solver->getAdvectionOperator()->getLocallyOwnedDofs());

        m_iterationNumber = m_solver->getIteration();
        m_time = m_solver->getTime();
        // get analytic and numeric values
        // TODO: only assign once (see. addAnalyticSolutionToOutput)
        m_solver->getAllAnalyticVelocities(m_solver->getTime(), ana_u,
                                           m_solver->m_supportPoints);
        m_solver->getAllAnalyticDensities(m_solver->getTime(), ana_rho,
                                          m_solver->m_supportPoints);

        //#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
        ana_rho.add(-1.0, num_rho);
        m_maxDensityError = ana_rho.linfty_norm();

        // calculate maximum analytic velocity norm
        const dealii::IndexSet& locally_owned_dofs =
                m_solver->getAdvectionOperator()->getLocallyOwnedDofs();
        m_maxUAnalytic = Math::maxVelocityNorm(ana_u,
                                               locally_owned_dofs);
        m_l2UAnalytic = Math::velocity2Norm(ana_u,
                                            locally_owned_dofs);

        // substract numeric from analytic velocity
        ana_u.at(0).add(-1.0, num_u.at(0));
        ana_u.at(1).add(-1.0, num_u.at(1));
        if (dim == 3) {
            ana_u.at(2).add(-1.0, num_u.at(2));
        }
        // calculate squares
        ana_u.at(0).scale(
                ana_u.at(0));
        ana_u.at(1).scale(
                ana_u.at(1));
        if (dim == 3) {
            ana_u.at(2).scale(
                    ana_u.at(2));
        }
        // calculate ||error (pointwise)||^2
        ana_u.at(0).add(ana_u.at(1));
        if (dim == 3) {
            ana_u.at(0).add(
                    ana_u.at(2));
        }

        // calculate || error (pointwise) ||
        //for all degrees of freedom on current processor
        dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
        dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
        for (; it != end; it++) {
            size_t i = *it;
            ana_u.at(0)(i) = sqrt(
                    ana_u.at(0)(i));
        }

        //rho
        const dealii::Function<dim> & f_rho =
                *m_solver->m_benchmark->getAnalyticRhoFunction(m_solver->getTime());
        dealii::Vector<double> local_errors(
                m_solver->getProblemDescription()->getMesh()->n_active_cells());
        dealii::VectorTools::integrate_difference(
                m_solver->getAdvectionOperator()->getMapping(),
                *m_solver->getAdvectionOperator()->getDoFHandler(), m_solver->getDensity(),
                f_rho, local_errors,
                *m_solver->getAdvectionOperator()->getQuadrature(),
                dealii::VectorTools::L2_norm);
        const double local_rho_error = local_errors.l2_norm();
        const double lsq = local_rho_error * local_rho_error;
        m_l2DensityError = sqrt(dealii::Utilities::MPI::sum(lsq,
                                                            MPI_COMM_WORLD));

        // u
        const dealii::Function<dim>& f_u =
                *m_solver->m_benchmark->getAnalyticUFunction(m_solver->getTime());
        AnalyticU<dim> ana_ux(f_u, 0);
        AnalyticU<dim> ana_uy(f_u, 1);
        // ux
        local_errors = 0;
        dealii::VectorTools::integrate_difference(
                m_solver->getAdvectionOperator()->getMapping(),
                *m_solver->getAdvectionOperator()->getDoFHandler(),
                m_solver->getVelocity().at(0), ana_ux, local_errors,
                *m_solver->getAdvectionOperator()->getQuadrature(),
                dealii::VectorTools::L2_norm);
        double one_component_local_error = local_errors.l2_norm();
        double total_local_error = one_component_local_error
                                   * one_component_local_error;
        // uy
        local_errors = 0;
        dealii::VectorTools::integrate_difference(
                m_solver->getAdvectionOperator()->getMapping(),
                *m_solver->getAdvectionOperator()->getDoFHandler(),
                m_solver->getVelocity().at(1), ana_uy, local_errors,
                *m_solver->getAdvectionOperator()->getQuadrature(),
                dealii::VectorTools::L2_norm);
        one_component_local_error = local_errors.l2_norm();
        total_local_error += one_component_local_error * one_component_local_error;
        // uz
        if (dim == 3) {
            AnalyticU<dim> ana_uz(f_u, 2);
            local_errors = 0;
            dealii::VectorTools::integrate_difference(
                    m_solver->getAdvectionOperator()->getMapping(),
                    *m_solver->getAdvectionOperator()->getDoFHandler(),
                    m_solver->getVelocity().at(2), ana_uz, local_errors,
                    *m_solver->getAdvectionOperator()->getQuadrature(),
                    dealii::VectorTools::L2_norm);
            one_component_local_error = local_errors.l2_norm();
            total_local_error += one_component_local_error
                                 * one_component_local_error;
        }
        const double total_global_error = sqrt(
                dealii::Utilities::MPI::sum(total_local_error,
                                            MPI_COMM_WORLD));

        m_maxVelocityError = ana_u.at(0).linfty_norm();
        m_l2VelocityError = total_global_error;

        CFDSolverUtilities::applyWriteableDensity(ana_rho,  m_solver->m_analyticDensity);
        CFDSolverUtilities::applyWriteableVelocity(ana_u,  m_solver->m_analyticVelocity);


    } /*update*/
    template void CompressibleErrorStats<2>::update();
    template void CompressibleErrorStats<3>::update();

    template<size_t dim>
    void CompressibleErrorStats<dim>::printHeaderLine() {
        if (is_MPI_rank_0()) {
            (*m_errorsTableFile)
                    << "#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
                    << endl;
        }
    }
    template void CompressibleErrorStats<2>::printHeaderLine();
    template void CompressibleErrorStats<3>::printHeaderLine();

    template<size_t dim>
    void CompressibleErrorStats<dim>::printNewLine() {
        if (not isUpToDate()) {
            update();
        }
        if (is_MPI_rank_0()) {
            (*m_errorsTableFile) << m_iterationNumber << " " << m_time << " "
                                 << m_maxUAnalytic << " " << m_maxVelocityError << " "
                                 << m_maxDensityError << " " << m_l2VelocityError << " "
                                 << m_l2DensityError << endl;
        }
    }
    template void CompressibleErrorStats<2>::printNewLine();
    template void CompressibleErrorStats<3>::printNewLine();

} /* namespace natrium */
