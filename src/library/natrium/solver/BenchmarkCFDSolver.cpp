/*
 * BenchmarkCFDSolver.cpp
 *
 *  Created on: 27.05.2014
 *      Author: kraemer
 */

#include "BenchmarkCFDSolver.h"

#include "deal.II/dofs/dof_tools.h"

namespace natrium {

template<size_t dim> BenchmarkCFDSolver<dim>::BenchmarkCFDSolver(
		boost::shared_ptr<SolverConfiguration> configuration,
		boost::shared_ptr<Benchmark<dim> > problemDescription) :
		CFDSolver<dim>(configuration, problemDescription), m_benchmark(
				problemDescription) {
	// initialize macroscopic variables
#ifdef WITH_TRILINOS_MPI
	m_analyticDensity.reinit(this->getDensity());
	for (size_t i = 0; i < dim; i++) {
		m_analyticVelocity.push_back(this->getVelocity().at(i));
	}
#else
	m_analyticDensity.reinit(this->getNumberOfDoFs(), true);
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
		m_errorStats = boost::make_shared<ErrorStats<dim> >(this, s.str());
	} else {
		m_errorStats = boost::make_shared<ErrorStats<dim> >(this);
	}

} /*BenchmarkCFDSolver constructor */
template BenchmarkCFDSolver<2>::BenchmarkCFDSolver(
		boost::shared_ptr<SolverConfiguration> configuration,
		boost::shared_ptr<Benchmark<2> > problemDescription);
template BenchmarkCFDSolver<3>::BenchmarkCFDSolver(
		boost::shared_ptr<SolverConfiguration> configuration,
		boost::shared_ptr<Benchmark<3> > problemDescription);

template<size_t dim>
void BenchmarkCFDSolver<dim>::output(size_t iteration) {
	CFDSolver<dim>::output(iteration);
	// output error table
	if (not this->getConfiguration()->isSwitchOutputOff()) {
		if (this->getIteration()
				% this->getConfiguration()->getOutputTableInterval() == 0) {
			m_errorStats->printNewLine();
		}
	}
} /*output*/
template void BenchmarkCFDSolver<2>::output(size_t iteration);
template void BenchmarkCFDSolver<3>::output(size_t iteration);

template<size_t dim>
void natrium::BenchmarkCFDSolver<dim>::getAllAnalyticDensities(double time,
		distributed_vector& analyticDensities,
		const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
	boost::shared_ptr<AdvectionOperator<dim> > adv_op = this->getAdvectionOperator();
	// get Function instance
	const boost::shared_ptr<dealii::Function<dim> > f_rho =
			m_benchmark->getAnalyticRhoFunction(time);
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
template
void natrium::BenchmarkCFDSolver<2>::getAllAnalyticDensities(double time,
		distributed_vector& analyticDensities,
		const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints) const;
template
void natrium::BenchmarkCFDSolver<3>::getAllAnalyticDensities(double time,
		distributed_vector& analyticDensities,
		const map<dealii::types::global_dof_index, dealii::Point<3> >& supportPoints) const;

template<size_t dim>
void natrium::BenchmarkCFDSolver<dim>::getAllAnalyticVelocities(double time,
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
template
void natrium::BenchmarkCFDSolver<2>::getAllAnalyticVelocities(double time,
		vector<distributed_vector>& analyticVelocities,
		const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints) const;
template
void natrium::BenchmarkCFDSolver<3>::getAllAnalyticVelocities(double time,
		vector<distributed_vector>& analyticVelocities,
		const map<dealii::types::global_dof_index, dealii::Point<3> >& supportPoints) const;

} /* namespace natrium */
