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
		shared_ptr<SolverConfiguration> configuration,
		shared_ptr<Benchmark<dim> > problemDescription) :
		CFDSolver<dim>(configuration, problemDescription), m_benchmark(
				problemDescription) {
	// initialize macroscopic variables
#ifdef WITH_TRILINOS_MPI
	m_supportPoints.resize(this->getAdvectionOperator()->getLocallyOwnedDofs().size());
	m_analyticDensity.reinit(
			this->getAdvectionOperator()->getLocallyOwnedDofs(),
			m_advectionOperator->getLocallyRelevantDofs(),
			MPI_COMM_WORLD);
	for (size_t i = 0; i < dim; i++) {
		m_analyticVelocity.push_back(
				distributed_vector(
						this->getAdvectionOperator()->getLocallyOwnedDofs(),
						m_advectionOperator->getLocallyRelevantDofs(),
						MPI_COMM_WORLD));
	}
#else
	m_supportPoints.resize(this->getNumberOfDoFs());
	m_analyticDensity.reinit(this->getNumberOfDoFs(), true);
	for (size_t i = 0; i < dim; i++) {
		m_analyticVelocity.push_back(
				distributed_vector(this->getNumberOfDoFs()));
	}
#endif
	dealii::DoFTools::map_dofs_to_support_points(
			this->getAdvectionOperator()->getMapping(),
			*(this->getAdvectionOperator()->getDoFHandler()), m_supportPoints);

// File for errors
	if ((not configuration->isSwitchOutputOff())) {
		std::stringstream s;
		s << configuration->getOutputDirectory().c_str() << "/errors_table.txt";
		m_errorStats = make_shared<ErrorStats<dim> >(this, s.str());
	} else {
		m_errorStats = make_shared<ErrorStats<dim> >(this);
	}

} /*BenchmarkCFDSolver constructor */
template BenchmarkCFDSolver<2>::BenchmarkCFDSolver(
		shared_ptr<SolverConfiguration> configuration,
		shared_ptr<Benchmark<2> > problemDescription);
template BenchmarkCFDSolver<3>::BenchmarkCFDSolver(
		shared_ptr<SolverConfiguration> configuration,
		shared_ptr<Benchmark<3> > problemDescription);

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

} /* namespace natrium */
