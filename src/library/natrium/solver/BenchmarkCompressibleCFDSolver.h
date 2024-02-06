//
// Created by dwilde3m on 25.02.20.
//

#ifndef NATRIUM_BENCHMARKCOMPRESSIBLECFDSOLVER_H
#define NATRIUM_BENCHMARKCOMPRESSIBLECFDSOLVER_H
#include "deal.II/numerics/data_out.h"

#include "CompressibleCFDSolver.h"
#include "SolverConfiguration.h"
#include "CompressibleErrorStats.h"
#include "BenchmarkCFDSolver.h"

#include "../problemdescription/CompressibleBenchmark.h"
#include "../utilities/BasicNames.h"
namespace natrium {
    template<size_t dim>
    class BenchmarkCompressibleCFDSolver : public CompressibleCFDSolver<dim>{
        template<size_t dim2> friend class CompressibleErrorStats;
    private:
        /// the problem description, pointed to explicitly as Benchmark
        boost::shared_ptr<CompressibleBenchmark<dim> > m_benchmark;

        /// support Point of DoFs, in the same order as the DoFs
        map<dealii::types::global_dof_index, dealii::Point<dim> > m_supportPoints;

        /// in order to save allocation time, allocate memory for analytic solution and errors
        distributed_vector m_analyticDensity;
        distributed_vector m_analyticTemperature;
        vector<distributed_vector> m_analyticVelocity;

        /// table out
        boost::shared_ptr<CompressibleErrorStats<dim> > m_compressibleErrorStats;

    public:
        /// constructor
        BenchmarkCompressibleCFDSolver(boost::shared_ptr<SolverConfiguration> configuration,
                           boost::shared_ptr<CompressibleBenchmark<dim> > problemDescription);

        /**
         * @short create output data and write to file
         */
        virtual void output(size_t iteration, bool is_final = false);

/// gives the possibility for Benchmark instances to add the analytic solution to output
        virtual void addAnalyticSolutionToOutput(dealii::DataOut<dim>& data_out) {
            distributed_vector rho;
            vector<distributed_vector> u;
            CFDSolverUtilities::getWriteableDensity(rho, m_analyticDensity,
                                                    this->getAdvectionOperator()->getLocallyOwnedDofs());
            CFDSolverUtilities::getWriteableVelocity(u, m_analyticVelocity,
                                                     this->getAdvectionOperator()->getLocallyOwnedDofs());
            getAllAnalyticVelocities(this->getTime(), u, m_supportPoints);
            getAllAnalyticDensities(this->getTime(), rho, m_supportPoints);
            CFDSolverUtilities::applyWriteableDensity(rho, m_analyticDensity);
            CFDSolverUtilities::applyWriteableVelocity(u, m_analyticVelocity);

            data_out.add_data_vector(m_analyticDensity, "rho_analytic");
            if (dim == 2) {
                data_out.add_data_vector(m_analyticVelocity.at(0), "ux_analytic");
                data_out.add_data_vector(m_analyticVelocity.at(1), "uy_analytic");
            } else { // dim == 3
                data_out.add_data_vector(m_analyticVelocity.at(0), "ux_analytic");
                data_out.add_data_vector(m_analyticVelocity.at(1), "uy_analytic");
                data_out.add_data_vector(m_analyticVelocity.at(2), "uz_analytic");
            }
        }

/// destructor
        virtual ~BenchmarkCompressibleCFDSolver() {

        }

        /**
         * @short get full analytic solution for the density field at time t
         */
        void getAllAnalyticDensities(double time,
                                     distributed_vector& analyticDensities,
                                     const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const;

        /**
         * @short get full analytic solution for the density field at time t
         */
        void getAllAnalyticVelocities(double time,
                                      vector<distributed_vector>& analyticVelocities,
                                      const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const;

        /**
         * @short set initial velocities
         * @param[out] initialVelocities vector of velocities; to be filled
         * @param[in] supportPoints the coordinates associated with each degree of freedom
         */
        void applyInitialVelocities(vector<distributed_vector>& initialVelocities,
                                    const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const;

//        const boost::shared_ptr<Benchmark<dim> > &getBenchmark() const {
//            const boost::shared_ptr<Benchmark<dim>> bench = ;
//            return m_benchmark;
//        }

        const map<dealii::types::global_dof_index, dealii::Point<dim> >& getSupportPoints() const {
            return m_supportPoints;
        }

        const boost::shared_ptr<CompressibleErrorStats<dim> >& getCompressibleErrorStats() const {
            return m_compressibleErrorStats;
        }
    }
    ;

} /* namespace natrium */

#endif //NATRIUM_BENCHMARKCOMPRESSIBLECFDSOLVER_H
