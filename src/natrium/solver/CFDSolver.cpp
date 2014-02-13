/**
 * @file CFDSolver.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CFDSolver.h"

#include <fstream>

#include "deal.II/numerics/data_out.h"

namespace natrium {

template<size_t dim>
CFDSolver<dim>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<dim> > problemDescription) {

	/// check if problem and solver configuration fit together
	configuration->checkProblem(problemDescription);
	m_problemDescription = problemDescription;
	m_configuration = configuration;

	/// Build boltzmann model
	if (Stencil_D2Q9 == configuration->getStencilType()) {
		m_boltzmannModel = make_shared<D2Q9IncompressibleModel>();
	}

	/// Build streaming data object
	if (Advection_SEDGMinLee == configuration->getAdvectionOperatorType()) {
		m_advectionOperator = make_shared<SEDGMinLee<dim> >(
				m_problemDescription->getTriangulation(),
				m_problemDescription->getBoundaries(),
				configuration->getOrderOfFiniteElement(), m_boltzmannModel,
				(configuration->getFluxType() == Flux_Central));
	}

	/// Calculate relaxation parameter and build collision model
	if (Collision_BGKTransformed == configuration->getCollisionType()) {
		double tau = BGKTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStep(), m_boltzmannModel);
		m_collisionModel = make_shared<BGKTransformed>(tau, m_boltzmannModel);
	}

	/// Build time integrator
	size_t numberOfDoFs = m_advectionOperator->getSystemMatrix().at(0).n();
	if (Integrator_RungeKutta5LowStorage
			== configuration->getTimeIntegratorType()) {
		m_timeIntegrator = make_shared<RungeKutta5LowStorage>(
				configuration->getTimeStep(), numberOfDoFs);
	}

	// initialize macroscopic variables
	vector<dealii::Point<dim> > supportPoints(numberOfDoFs);
	m_advectionOperator->mapDoFsToSupportPoints(supportPoints);
	m_density.reinit(numberOfDoFs);
	for (size_t i = 0; i < dim; i++) {
		distributed_vector vi(numberOfDoFs);
		m_velocity.push_back(vi);
	}
	m_problemDescription->applyInitialDensities(m_density, supportPoints);
	m_problemDescription->applyInitialVelocities(m_velocity, supportPoints);

	// Initialize distribution functions
	for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
		distributed_vector fi(m_advectionOperator->getNumberOfDoFs());
		m_f.push_back(fi);
	}
	vector<double> feq(m_boltzmannModel->getQ());
	numeric_vector u(dim);
	for (size_t i = 0; i < numberOfDoFs; i++) {
		for (size_t j = 0; j < dim; j++) {
			u(j) = m_velocity.at(j)(i);
		}
		m_boltzmannModel->getEquilibriumDistributions(feq, u, m_density(i));
		for (size_t j = 0; j < m_boltzmannModel->getQ(); j++) {
			m_f.at(j)(i) = feq.at(j);
		}
	}

} /* Constructor */
template CFDSolver<2>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<2> > problemDescription);
template CFDSolver<3>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<3> > problemDescription);

template<size_t dim>
void CFDSolver<dim>::stream() {
	const vector<distributed_sparse_matrix>& systemMatrix =
			m_advectionOperator->getSystemMatrix();
	// no streaming in direction 0; begin with 1
	for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
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
		if (i % 100 == 0) {
			cout << "Iteration " << i << endl;
		}
		stream();
		collide();
		output(i);
	}
}
template void CFDSolver<2>::run();
template void CFDSolver<3>::run();

template<size_t dim>
void CFDSolver<dim>::output(size_t iteration) {
	std::stringstream str;
	str << m_configuration->getOutputDirectory().c_str() << "/t_" << iteration
			<< ".vtu";
	std::string filename = str.str();
	std::ofstream vtu_output(filename.c_str());
	dealii::DataOut<dim> data_out;
	data_out.attach_dof_handler(*m_advectionOperator->getDoFHandler());
	data_out.add_data_vector(m_density, "rho");
	for (size_t i = 0; i < dim; i++) {
		std::stringstream vi;
		vi << "v_" << i;
		data_out.add_data_vector(m_velocity.at(i), vi.str().c_str());
	}
	data_out.build_patches();
	data_out.write_vtu(vtu_output);
}
template void CFDSolver<2>::output(size_t iteration);
template void CFDSolver<3>::output(size_t iteration);

} /* namespace natrium */
