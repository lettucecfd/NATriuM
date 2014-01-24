/**
 * @file CFDSolver.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CFDSolver.h"

namespace natrium {

template<size_t dim>
CFDSolver<dim>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<dim> > problemDescription) {

	/// check if problem and solver configuration fit together
	configuration->checkProblem(problemDescription);
	m_problemDescription = problemDescription;
	m_configuration = configuration;

	/// Build streaming data object
	if (Advection_SEDGMinLee == configuration->getAdvectionOperatorType()) {
		/*m_streamingData = make_shared<DataMinLee2011<dim> >(
		 m_problemDescription->getTriangulation(),
		 configuration->getOrderOfFiniteElement());
		 */}

	/// Build boltzmann model
	if (Stencil_D2Q9 == configuration->getStencilType()) {
		m_boltzmannModel = make_shared<D2Q9IncompressibleModel>();
	}

	/// Build collision model
	if (Collision_BGKTransformed == configuration->getCollisionType()) {
		m_collisionModel = make_shared<BGKTransformed>(
				m_problemDescription->getRelaxationParameter(),
				m_boltzmannModel);
	}

	/// Build time integrator
	/*if (Integrator_RungeKutta5LowStorage == configuration->getTimeIntegratorType()){
	 m_timeIntegrator = make_shared<RungeKutta5LowStorage>();
	 }*/

} /* Constructor */
template CFDSolver<2>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<2> > problemDescription);
template CFDSolver<3>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<3> > problemDescription);

template<size_t dim>
void CFDSolver<dim>::stream() {
	const vector<distributed_sparse_matrix>& systemMatrix =
			m_advectionOperator->getSystemMatrix();
	for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
		m_timeIntegrator->step(m_f.at(i), systemMatrix.at(i));
	}

}
template void CFDSolver<2>::stream();
template void CFDSolver<3>::stream();

template<size_t dim>
void CFDSolver<dim>::collide() {
	m_collisionModel->collideAll(m_f, m_density, m_velocity);
}
template void CFDSolver<2>::collide();
template void CFDSolver<3>::collide();

template<size_t dim>
void CFDSolver<dim>::reassemble() {
	m_advectionOperator->reassemble();
}
template void CFDSolver<2>::reassemble();
template void CFDSolver<3>::reassemble();

template<size_t dim>
void CFDSolver<dim>::run() {
	size_t N = m_configuration->getNumberOfTimeSteps();
	for (size_t i = 0; i < N; i++) {
		stream();
		collide();
	}
}
template void CFDSolver<2>::run();
template void CFDSolver<3>::run();

} /* namespace natrium */
