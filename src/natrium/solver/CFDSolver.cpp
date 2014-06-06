/**
 * @file CFDSolver.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CFDSolver.h"

#include <fstream>
#include <sstream>

#include <boost/filesystem.hpp>

#include "deal.II/numerics/data_out.h"
#include "deal.II/fe/component_mask.h"
#include "deal.II/base/logstream.h"
#include "deal.II/grid/grid_tools.h"

#include "PhysicalProperties.h"

#include "../problemdescription/BoundaryCollection.h"

#include "../utilities/Logging.h"

namespace natrium {

template<size_t dim>
CFDSolver<dim>::CFDSolver(shared_ptr<SolverConfiguration> configuration,
		shared_ptr<ProblemDescription<dim> > problemDescription) {

	/// Create output directory
	if (not configuration->isSwitchOutputOff()) {
		configuration->prepareOutputDirectory();
	}

	// Create file for errors
	if ((not configuration->isSwitchOutputOff())
	/*and (configuration->getOutputTableInterval()
	 < configuration->getNumberOfTimeSteps())*/) {
		std::stringstream s;
		s << configuration->getOutputDirectory().c_str()
				<< "/results_table.txt";
		if (this->getIterationStart() > 0) {
			m_tableFile = make_shared<std::fstream>(s.str().c_str(),
					std::fstream::out | std::fstream::app);
		} else {
			m_tableFile = make_shared<std::fstream>(s.str().c_str(),
					std::fstream::out);
			(*m_tableFile) << "#  i      t      max |u_numeric|    kinE"
					<< endl;
		}
	}

	// CONFIGURE LOGGER
	if (configuration->isSwitchOutputOff()) {
		LOGGER().setConsoleLevel(SILENT);
		LOGGER().setFileLevel(SILENT);
	} else {
		LOGGER().setConsoleLevel(
				LogLevel(configuration->getCommandLineVerbosity()));
		LOGGER().setFileLevel(ALL);
		LOGGER().setLogFile(
				(boost::filesystem::path(configuration->getOutputDirectory())
						/ "natrium.log").string());
	}

	/// check if problem and solver configuration fit together
	configuration->checkProblem(problemDescription);
	m_problemDescription = problemDescription;
	m_configuration = configuration;

	/// Build boltzmann model
	if (Stencil_D2Q9 == configuration->getStencil()) {
		m_boltzmannModel = make_shared<D2Q9IncompressibleModel>(
				configuration->getStencilScaling());
	}

	/// Build streaming data object
	if (SEDG == configuration->getAdvectionScheme()) {
		// if this string is "": matrix is reassembled and not read from file
		string whereAreTheStoredMatrices;
		if (configuration->isRestartAtLastCheckpoint()) {
			whereAreTheStoredMatrices = configuration->getOutputDirectory();
		}
		// create SEDG MinLee by reading the system matrices from files or assembling
		// TODO estimate time for assembly
		m_advectionOperator = make_shared<SEDGMinLee<dim> >(
				m_problemDescription->getTriangulation(),
				m_problemDescription->getBoundaries(),
				configuration->getSedgOrderOfFiniteElement(), m_boltzmannModel,
				whereAreTheStoredMatrices,
				(CENTRAL == configuration->getSedgFluxType()));
	}

	if (configuration->isRestartAtLastCheckpoint()) {
		// read iteration number and time from file
		std::stringstream filename;
		filename << m_configuration->getOutputDirectory() << "/checkpoint.dat";
		std::ifstream ifile(filename.str().c_str());
		ifile >> m_iterationStart;
		ifile >> m_time;
		// print out message
		LOG(BASIC) << "Restart at iteration " << m_iterationStart
				<< " (simulation time = " << m_time << " s)." << endl;
	} else {
		m_iterationStart = 0;
		m_time = 0;
	}

/// Calculate relaxation parameter and build collision model
	double tau = 0.0;
	if (BGK_WITH_TRANSFORMED_DISTRIBUTION_FUNCTIONS
			== configuration->getCollisionScheme()) {
		tau = BGKTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), m_boltzmannModel);
		m_collisionModel = make_shared<BGKTransformed>(tau, m_boltzmannModel);
	}

/// Build time integrator
	size_t numberOfDoFs = m_advectionOperator->getNumberOfDoFs();
	if (RUNGE_KUTTA_5STAGE == configuration->getTimeIntegrator()) {
		m_timeIntegrator = make_shared<
				RungeKutta5LowStorage<distributed_sparse_block_matrix,
						distributed_block_vector> >(
				configuration->getTimeStepSize(), numberOfDoFs,
				m_boltzmannModel->getQ() - 1);
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
	double charU = problemDescription->getCharacteristicVelocity();
	if (charU == 0.0) {
		charU = maxU;
	}
	double dx = dealii::GridTools::minimal_cell_diameter(
			*(m_problemDescription->getTriangulation()));
	LOG(WELCOME) << "------ NATriuM solver ------" << endl;
	LOG(WELCOME) << "viscosity:       " << problemDescription->getViscosity()
			<< " m^2/s" << endl;
	LOG(WELCOME) << "char. length:    "
			<< problemDescription->getCharacteristicLength() << " m" << endl;
	LOG(WELCOME) << "max |u_0|:       "
			<< maxU * problemDescription->getCharacteristicLength() << " m/s"
			<< endl;
	LOG(WELCOME) << "Reynolds number: "
			<< (charU * problemDescription->getCharacteristicLength())
					/ problemDescription->getViscosity() << endl;
	LOG(WELCOME) << "Mach number:     "
				<< charU / m_boltzmannModel->getSpeedOfSound() << endl;
	LOG(WELCOME) << "Stencil scaling: "
				<< configuration->getStencilScaling() << endl;
	LOG(WELCOME) << "Recommended dt:  "
			<< m_collisionModel->calculateOptimalTimeStep(dx, m_boltzmannModel)
			<< " s" << endl;
	LOG(WELCOME) << "Actual dt:       " << configuration->getTimeStepSize()
			<< " s" << endl;
	LOG(WELCOME) << "dx:              " << dx << endl;
	LOG(WELCOME) << "tau:             " << tau << endl;
	LOG(WELCOME) << "----------------------------" << endl;

// Initialize distribution functions
	if (configuration->isRestartAtLastCheckpoint()) {
		loadDistributionFunctionsFromFiles(
				m_configuration->getOutputDirectory());
	} else {
		initializeDistributions();
	}

	// initialize boundary dof indicator
	std::set<dealii::types::boundary_id> boundaryIndicators;
	typename BoundaryCollection<dim>::ConstIterator it =
			m_problemDescription->getBoundaries()->getBoundaries().begin();
	for (; it != m_problemDescription->getBoundaries()->getBoundaries().end();
			it++) {
		if (not it->second->isPeriodic()) {
			boundaryIndicators.insert(it->first);
		}
	}
	m_isDoFAtBoundary.resize(getNumberOfDoFs());
	dealii::DoFTools::extract_dofs_with_support_on_boundary(
			*(m_advectionOperator->getDoFHandler()), dealii::ComponentMask(),
			m_isDoFAtBoundary, boundaryIndicators);
	size_t nofBoundaryNodes = 0;
	for (size_t i = 0; i < m_isDoFAtBoundary.size(); i++) {
		if (m_isDoFAtBoundary.at(i)) {
			nofBoundaryNodes += 1;
		}
	}
	LOG(DETAILED) << "Number of non-periodic boundary dofs: 9*"
			<< nofBoundaryNodes << endl;

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
	//f.print(cout);
	f.add(m_timeIntegrator->getTimeStepSize(),
			m_advectionOperator->getSystemVector());
	m_time += m_timeIntegrator->getTimeStepSize();
	//f.add(1.0, m_advectionOperator->getSystemVector());
	//m_advectionOperator->getSystemVector().print(cout);
	//f.print(cout);
	//assert (false);
	/*
	 for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
	 m_timeIntegrator->step(m_f.at(i), systemMatrix.block(i-1,i-1));
	 }*/

}
template void CFDSolver<2>::stream();
template void CFDSolver<3>::stream();

template<size_t dim>
void CFDSolver<dim>::collide() {
	//m_collisionModel->collideAll(m_f, m_density, m_velocity, m_isBoundary);
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
	// TODO estimate runtime
	for (m_i = m_iterationStart; m_i <= N; m_i++) {
		output(m_i);
		stream();
		collide();
	}
}
template void CFDSolver<2>::run();
template void CFDSolver<3>::run();

template<size_t dim>
void CFDSolver<dim>::output(size_t iteration) {
	// output: vector fields as .vtu files
	if (not m_configuration->isSwitchOutputOff()) {
		if (iteration % 100 == 0) {
			LOG(BASIC) << "Iteration " << iteration << ",  t = " << m_time
					<< endl;
		}
		if (iteration % m_configuration->getOutputSolutionInterval() == 0) {
			std::stringstream str;
			str << m_configuration->getOutputDirectory().c_str() << "/t_"
					<< iteration << ".vtu";
			std::string filename = str.str();
			std::ofstream vtu_output(filename.c_str());
			dealii::DataOut<dim> data_out;
			data_out.attach_dof_handler(*m_advectionOperator->getDoFHandler());
			data_out.add_data_vector(m_density, "rho");
			if (dim == 2) {
				data_out.add_data_vector(m_velocity.at(0), "ux");
				data_out.add_data_vector(m_velocity.at(1), "uy");
			} else { //dim == 3
				data_out.add_data_vector(m_velocity.at(0), "ux");
				data_out.add_data_vector(m_velocity.at(1), "uy");
				data_out.add_data_vector(m_velocity.at(2), "uz");
			}
			addAnalyticSolutionToOutput(data_out);
			/// For Benchmarks: add analytic solution
			data_out.build_patches();
			data_out.write_vtu(vtu_output);
		}

		// output: table
		if (iteration % m_configuration->getOutputTableInterval() == 0) {
			//#  i      t         max |u_numeric|  kinE  mass"
			(*m_tableFile) << getIteration() << " " << getTime()
					<< " " << getMaxVelocityNorm() << " "
					<< PhysicalProperties<dim>::kineticEnergy(getVelocity(),
							getDensity()) << endl;
			//<< " " << PhysicalProperties<dim>::mass(getVelocity(), getDensity());
		}

		// output: checkpoint
		if (iteration % m_configuration->getOutputCheckpointInterval() == 0) {
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
			// time
			outfile << m_time << endl;
		}
	}
}
template void CFDSolver<2>::output(size_t iteration);
template void CFDSolver<3>::output(size_t iteration);

template<size_t dim>
void CFDSolver<dim>::initializeDistributions() {
	LOG(BASIC) << "Initialize distribution functions: ";
	vector<double> feq(m_boltzmannModel->getQ());
	numeric_vector u(dim);

	// Initialize f with the equilibrium distribution functions
	m_f.reinit(m_boltzmannModel->getQ(),
			m_advectionOperator->getNumberOfDoFs());
	for (size_t i = 0; i < m_velocity.at(0).size(); i++) {
		for (size_t j = 0; j < dim; j++) {
			u(j) = m_velocity.at(j)(i);
		}
		m_boltzmannModel->getEquilibriumDistributions(feq, u, m_density(i));
		for (size_t j = 0; j < m_boltzmannModel->getQ(); j++) {
			m_f.at(j)(i) = feq.at(j);
		}
	}

	switch (m_configuration->getInitializationScheme()) {
	case EQUILIBRIUM: {
		LOG(BASIC) << "Equilibrium distribution functions" << endl;
		// do nothing else
		break;
	}
	case ITERATIVE: {
		LOG(BASIC) << "Iterative procedure" << endl;
		LOG(DETAILED) << "residual = "
				<< m_configuration->getIterativeInitializationResidual();
		LOG(DETAILED) << ", max iterations = "
				<< m_configuration->getIterativeInitializationNumberOfIterations()
				<< endl;
		// Iterative procedure; leading to consistent initial values
		size_t loopCount = 0;
		double residual = 10;
		const bool inInitializationProcedure = true;
		distributed_vector oldDensities;
		while (residual > m_configuration->getIterativeInitializationResidual()) {
			if (loopCount
					> m_configuration->getIterativeInitializationNumberOfIterations()) {
				throw CFDSolverException(
						"The iterative Initialization of equilibrium distribution functions failed. Either use another initialization procedure or change the preferences for iterative initialization (residual or maximal number of init iterations).");
			}
			oldDensities = m_density;
			stream();
			// collide without recalculating velocities
			m_collisionModel->collideAll(m_f, m_density, m_velocity,
					m_isDoFAtBoundary, inInitializationProcedure);
			oldDensities -= m_density;
			residual = oldDensities.norm_sqr();
			loopCount++;
		}
		LOG(DETAILED) << "Residual " << residual << " reached after "
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

	LOG(BASIC) << "Initialize distribution functions: done." << endl;
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
	m_f.reinit(m_boltzmannModel->getQ(),
			m_advectionOperator->getNumberOfDoFs());
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

