//
// Created by dwilde3m on 25.02.20.
//

#include "BenchmarkCompressibleCFDSolver.h"
#include "deal.II/dofs/dof_tools.h"
namespace natrium {
    template<size_t dim>
    BenchmarkCompressibleCFDSolver<dim>::BenchmarkCompressibleCFDSolver(
           boost::shared_ptr<SolverConfiguration> configuration,
            boost::shared_ptr<CompressibleBenchmark<dim> > problemDescription) :
            CompressibleCFDSolver<dim>(configuration, problemDescription), m_benchmark(
            problemDescription) {
        // initialize macroscopic variables
#ifdef WITH_TRILINOS_MPI
        m_analyticDensity.reinit(this->getDensity());
        m_analyticTemperature.reinit(this->getTemperature());
        for (size_t i = 0; i < dim; i++) {
            m_analyticVelocity.push_back(this->getVelocity().at(i));
        }
#else
        m_analyticDensity.reinit(this->getNumberOfDoFs(), true);
        m_analyticTemperature.reinit(this->getNumberOfDoFs(), true);
	for (size_t i = 0; i < dim; i++) {
		m_analyticVelocity.push_back(
				distributed_vector(this->getNumberOfDoFs()));
	}
#endif
        this->getAdvectionOperator()->mapDoFsToSupportPoints(m_supportPoints);

// File for errors
        if ((not configuration->isSwitchOutputOff())) {
            std::stringstream s;
            s << configuration->getOutputDirectory().c_str() << "/errors_table.txt";
            m_compressibleErrorStats = boost::make_shared<CompressibleErrorStats<dim> >(this, s.str());
        } else {
            m_compressibleErrorStats = boost::make_shared<CompressibleErrorStats<dim> >(this);
        }

    } /*BenchmarkCFDSolver constructor */

    template<size_t dim>
    void BenchmarkCompressibleCFDSolver<dim>::output(size_t iteration, bool ) {
        CompressibleCFDSolver<dim>::output(iteration);
        // output error table
        if (not this->getConfiguration()->isSwitchOutputOff()) {
            if (this->getIteration()
                % this->getConfiguration()->getOutputTableInterval() == 0) {
                m_compressibleErrorStats->printNewLine();
            }
        }
    } /*output*/

    template<size_t dim>
    void BenchmarkCompressibleCFDSolver<dim>::getAllAnalyticDensities(double time,
                                                          distributed_vector& analyticDensities,
                                                          const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
        boost::shared_ptr<AdvectionOperator<dim> > adv_op = this->getAdvectionOperator();
        // get Function instance
        const boost::shared_ptr<dealii::Function<dim> > f_rho =
                m_benchmark->getAnalyticRhoFunction(time);
        cout << time << endl;
        const unsigned int dofs_per_cell = adv_op->getFe()->dofs_per_cell;
        vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
        typename dealii::DoFHandler<dim>::active_cell_iterator cell =
                adv_op->getDoFHandler()->begin_active(), endc =
                adv_op->getDoFHandler()->end();
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                cell->get_dof_indices(local_dof_indices);
                for (size_t i = 0; i < dofs_per_cell; i++) {
                    assert(
                            supportPoints.find(local_dof_indices.at(i))
                            != supportPoints.end());
                    analyticDensities(local_dof_indices.at(i)) = f_rho->value(
                            supportPoints.at(local_dof_indices.at(i)));
                }
            } /* if is locally owned */
        } /* for all cells */
    }

    template<size_t dim>
    void BenchmarkCompressibleCFDSolver<dim>::getAllAnalyticVelocities(double time,
                                                           vector<distributed_vector>& analyticVelocities,
                                                           const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
        boost::shared_ptr<AdvectionOperator<dim> > adv_op = this->getAdvectionOperator();
        // get Function instance
        const boost::shared_ptr<dealii::Function<dim> > f_u =
                m_benchmark->getAnalyticUFunction(time);
        const unsigned int dofs_per_cell = adv_op->getFe()->dofs_per_cell;
        vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
        typename dealii::DoFHandler<dim>::active_cell_iterator cell =
                adv_op->getDoFHandler()->begin_active(), endc =
                adv_op->getDoFHandler()->end();
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                cell->get_dof_indices(local_dof_indices);
                for (size_t i = 0; i < dofs_per_cell; i++) {
                    assert(
                            supportPoints.find(local_dof_indices.at(i))
                            != supportPoints.end());
                    for (size_t component = 0; component < dim; component++) {
                        analyticVelocities.at(component)(local_dof_indices.at(i)) =
                                f_u->value(
                                        supportPoints.at(local_dof_indices.at(i)),
                                        component);
                    }
                }
            } /* if is locally owned */
        } /* for all cells */
    }


// Explicit instatiation
    template class BenchmarkCompressibleCFDSolver<2>;
    template class BenchmarkCompressibleCFDSolver<3>;

} /* namespace natrium */