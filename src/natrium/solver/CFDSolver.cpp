/**
 * @file CFDSolver.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CFDSolver.h"

#include <fstream>
#include <sstream>

#include "deal.II/numerics/data_out.h"

#include "../utilities/Logging.h"

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
		m_boltzmannModel = make_shared<D2Q9IncompressibleModel>(
				configuration->getDQScaling());
	}

	/// Build streaming data object
	if (Advection_SEDGMinLee == configuration->getAdvectionOperatorType()) {
		// if this string is "": matrix is reassembled and not read from file
		string whereAreTheStoredMatrices;
		if (configuration->isRestart()) {
			whereAreTheStoredMatrices = configuration->getOutputDirectory();
		}
		// create SEDG MinLee by reading the system matrices from files or assembling
		m_advectionOperator = make_shared<SEDGMinLee<dim> >(
				m_problemDescription->getTriangulation(),
				m_problemDescription->getBoundaries(),
				configuration->getOrderOfFiniteElement(), m_boltzmannModel,
				whereAreTheStoredMatrices,
				(configuration->getFluxType() == Flux_Central));
	}

	if (configuration->isRestart()) {
		// read iteration number from file
		std::stringstream filename;
		filename << m_configuration->getOutputDirectory() << "/checkpoint.dat";
		std::ifstream ifile(filename.str().c_str());
		ifile >> m_iterationStart;
		// print out message
		*(Logging::BASIC) << "Restart at iteration " << m_iterationStart
				<< endl;
	} else {
		m_iterationStart = 0;
	}

/// Calculate relaxation parameter and build collision model
	double tau = 0.0;
	if (Collision_BGKTransformed == configuration->getCollisionType()) {
		tau = BGKTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStep(), m_boltzmannModel);
		m_collisionModel = make_shared<BGKTransformed>(tau, m_boltzmannModel);
	}

/// Build time integrator
	size_t numberOfDoFs = m_advectionOperator->getNumberOfDoFs();
	if (Integrator_RungeKutta5LowStorage
			== configuration->getTimeIntegratorType()) {
		m_timeIntegrator = make_shared<RungeKutta5LowStorage<distributed_sparse_block_matrix, distributed_block_vector> >(
				configuration->getTimeStep(), numberOfDoFs, m_boltzmannModel->getQ() - 1);
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

// OUTPUT
	double maxU = getMaxVelocityNorm();
	*(Logging::BASIC) << "------ NATriuM solver ------" << endl;
	*(Logging::BASIC) << "viscosity:       "
			<< problemDescription->getViscosity() << " m^2/s" << endl;
	*(Logging::BASIC) << "char. length:    "
			<< problemDescription->getCharacteristicLength() << " m" << endl;
	*(Logging::BASIC) << "max |u_0|:       "
			<< maxU * problemDescription->getCharacteristicLength() << " m/s"
			<< endl;
	*(Logging::BASIC) << "Reynolds number: "
			<< (maxU * problemDescription->getCharacteristicLength())
					/ problemDescription->getViscosity() << endl;
	*(Logging::BASIC) << "Recommended dt:  "
			<< m_collisionModel->calculateOptimalTimeStep(
					problemDescription->getViscosity(), m_boltzmannModel)
			<< " s" << endl;
	*(Logging::BASIC) << "Actual dt:       " << configuration->getTimeStep()
			<< " s" << endl;
	*(Logging::BASIC) << "tau:             " << tau << endl;
	*(Logging::BASIC) << "----------------------------" << endl;

// Initialize distribution functions
	if (configuration->isRestart()) {
		loadDistributionFunctionsFromFiles(
				m_configuration->getOutputDirectory());
	} else {
		initializeDistributions();
	}

}
/* Constructor */
template CFDSolver<2>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<2> > problemDescription);
template CFDSolver<3>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<3> > problemDescription);

template<size_t dim>
void CFDSolver<dim>::stream() {
	const distributed_sparse_block_matrix& systemMatrix =
			m_advectionOperator->getSystemMatrix();
	distributed_block_vector& f = m_f.getFStream();
// no streaming in direction 0; begin with 1
	m_timeIntegrator->step(f, systemMatrix);
	f.add(m_timeIntegrator->getTimeStepSize(), m_advectionOperator->getSystemVector());
	/*
	for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
		m_timeIntegrator->step(m_f.at(i), systemMatrix.block(i-1,i-1));
	}*/

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
	for (size_t i = m_iterationStart; i < N; i++) {
		if (i % 100 == 0) {
			*(Logging::BASIC) << "Iteration " << i << endl;
		}
		output(i);
		stream();
		collide();
	}
}
template void CFDSolver<2>::run();
template void CFDSolver<3>::run();

template<size_t dim>
void CFDSolver<dim>::output(size_t iteration) {
	// output: vector fields as .vtu files
	if ((out_VectorFields & m_configuration->getOutputFlags()) != 0) {
		if (iteration % m_configuration->getOutputVectorFieldsEvery() == 0) {
			std::stringstream str;
			str << m_configuration->getOutputDirectory().c_str() << "/t_"
					<< iteration << ".vtu";
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
	}
	// output: checkpoint
	if ((out_Checkpoints & m_configuration->getOutputFlags()) != 0) {
		if (iteration % m_configuration->getOutputCheckpointEvery() == 0) {
			// advection matrices
			m_advectionOperator->saveCheckpoint(
					m_configuration->getOutputDirectory());
			// distribution functions
			saveDistributionFunctionsToFiles(
					m_configuration->getOutputDirectory());
			// iteration
			std::stringstream filename;
			filename << m_configuration->getOutputDirectory()
					<< "/checkpoint.dat";
			std::ofstream outfile(filename.str().c_str());
			outfile << iteration << endl;
		}
	}
}
template void CFDSolver<2>::output(size_t iteration);
template void CFDSolver<3>::output(size_t iteration);

template<size_t dim>
void CFDSolver<dim>::initializeDistributions() {
	(*Logging::BASIC) << "Initialize distribution functions: ";
	vector<double> feq(m_boltzmannModel->getQ());
	numeric_vector u(dim);

	// Initialize f with the equilibrium distribution functions
	m_f.reinit(m_boltzmannModel->getQ(), m_advectionOperator->getNumberOfDoFs());
	for (size_t i = 0; i < m_velocity.at(0).size(); i++) {
		for (size_t j = 0; j < dim; j++) {
			u(j) = m_velocity.at(j)(i);
		}
		m_boltzmannModel->getEquilibriumDistributions(feq, u, m_density(i));
		for (size_t j = 0; j < m_boltzmannModel->getQ(); j++) {
			m_f.at(j)(i) = feq.at(j);
		}
	}

	switch (m_configuration->getDistributionInitType()) {
	case Equilibrium: {
		(*Logging::BASIC) << "Equilibrium distribution functions" << endl;
		// do nothing else
		break;
	}
	case Iterative: {
		(*Logging::BASIC) << "Iterative procedure" << endl;
		(*Logging::FULL) << "residual = "
				<< m_configuration->getStopDistributionInitResidual();
		(*Logging::FULL) << ", max iterations = "
				<< m_configuration->getMaxDistributionInitIterations() << endl;
		// Iterative procedure; leading to consistent initial values
		size_t loopCount = 0;
		double residual = 10;
		const bool inInitializationProcedure = true;
		distributed_vector oldDensities;
		while (residual > m_configuration->getStopDistributionInitResidual()) {
			if (loopCount
					> m_configuration->getMaxDistributionInitIterations()) {
				throw CFDSolverException(
						"The iterative Initialization of equilibrium distribution functions failed. Either use another initialization procedure or change the preferences for iterative initialization (residual or maximal number of init iterations).");
			}
			oldDensities = m_density;
			stream();
			// collide without recalculating velocities
			m_collisionModel->collideAll(m_f, m_density, m_velocity,
					inInitializationProcedure);
			oldDensities -= m_density;
			residual = oldDensities.norm_sqr();
			loopCount++;
		}
		(*Logging::FULL) << "Residual " << residual << " reached after "
				<< loopCount << " iterations." << endl;

		for (size_t i = 0; i < m_velocity.at(0).size(); i++) {
			for (size_t j = 0; j < dim; j++) {
				u(j) = m_velocity.at(j)(i);
			}
			m_boltzmannModel->getEquilibriumDistributions(feq, u, m_density(i));
			for (size_t j = 0; j < m_boltzmannModel->getQ(); j++) {
				m_f.at(j)(i) = feq.at(j);
			}
		}
		break;
	}
	default: {
		throw CFDSolverException(
				"Error in CFDSolver::InitializeDistributions. A part of the code was reached, which should never be reached.");
		break;
	}
	}

	(*Logging::BASIC) << "Initialize distribution functions: done." << endl;
}

template void CFDSolver<2>::initializeDistributions();
template void CFDSolver<3>::initializeDistributions();

template<size_t dim>
void CFDSolver<dim>::saveDistributionFunctionsToFiles(const string& directory) {
	for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
		// filename
		std::stringstream filename;
		filename << directory << "/checkpoint_f_" << i << ".dat";
		std::ofstream file(filename.str().c_str());
		m_f.at(i).block_write(file);
	}
}
template void CFDSolver<2>::saveDistributionFunctionsToFiles(
		const string& directory);
template void CFDSolver<3>::saveDistributionFunctionsToFiles(
		const string& directory);

template<size_t dim>
void CFDSolver<dim>::loadDistributionFunctionsFromFiles(
		const string& directory) {
	// create vectors
	m_f.reinit(m_boltzmannModel->getQ(), m_advectionOperator->getNumberOfDoFs());
	// read the distribution functions from file
	try {
		for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
			// filename
			std::stringstream filename;
			filename << directory << "/checkpoint_f_" << i << ".dat";
			std::ifstream file(filename.str().c_str());
			m_f.at(i).block_read(file);
		}
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw CFDSolverException(
				"An error occurred while reading the distribution functions from file: Please switch off the restart option to start the simulation from the beginning.");
	}
}
template void CFDSolver<2>::loadDistributionFunctionsFromFiles(
		const string& directory);
template void CFDSolver<3>::loadDistributionFunctionsFromFiles(
		const string& directory);

} /* namespace natrium */

