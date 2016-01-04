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
#include "deal.II/lac/constraint_matrix.h"
#include "deal.II/base/index_set.h"

#include "PhysicalProperties.h"

#include "../stencils/D2Q9.h"
#include "../stencils/D3Q19.h"
#include "../stencils/D3Q15.h"
#include "../stencils/D3Q27.h"
#include "../stencils/Stencil.h"

#include "../collision/BGKStandard.h"
#include "../collision/BGKStandardTransformed.h"
#include "../collision/BGKSteadyState.h"
#include "../collision/BGKIncompressible.h"
#include "../collision/MRTStandard.h"
#include "../collision/KBCStandard.h"

#include "../problemdescription/BoundaryCollection.h"
#include "../problemdescription/NonlinearBoundary.h"

#include "../timeintegration/ThetaMethod.h"
#include "../timeintegration/RungeKutta5LowStorage.h"
#include "../timeintegration/ExponentialTimeIntegrator.h"
#include "../timeintegration/DealIIWrapper.h"

#include "../utilities/Logging.h"
#include "../utilities/CFDSolverUtilities.h"
#include "../utilities/MPIGuard.h"
#include "../utilities/Info.h"

namespace natrium {

template<size_t dim>
CFDSolver<dim>::CFDSolver(boost::shared_ptr<SolverConfiguration> configuration,
		boost::shared_ptr<ProblemDescription<dim> > problemDescription) {

	/// Create output directory
	if (not configuration->isSwitchOutputOff()) {
		try {
			configuration->prepareOutputDirectory();
		} catch (ConfigurationException & e) {
			natrium_errorexit(e.what());
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

	/// check if problem's boundary conditions are well defined
	bool boundaries_ok = problemDescription->checkBoundaryConditions();
	if (!boundaries_ok) {
		natrium_errorexit("Boundary conditions do no fit to triangulation.");
	}

	/// check if problem and solver configuration fit together
	try {
		configuration->checkProblem(problemDescription);
	} catch (ConfigurationException & e) {
		natrium_errorexit(e.what());
	}
	m_problemDescription = problemDescription;
	m_configuration = configuration;

	/// Build boltzmann model
	if (Stencil_D2Q9 == configuration->getStencil()) {
		m_stencil = boost::make_shared<D2Q9>(
				configuration->getStencilScaling());
	} else if (Stencil_D3Q19 == configuration->getStencil()) {
		m_stencil = boost::make_shared<D3Q19>(
				configuration->getStencilScaling());
	} else if (Stencil_D3Q15 == configuration->getStencil()) {
		m_stencil = boost::make_shared<D3Q15>(
				configuration->getStencilScaling());
	} else if (Stencil_D3Q27 == configuration->getStencil()) {
		m_stencil = boost::make_shared<D3Q27>(
				configuration->getStencilScaling());
	} else {
		natrium_errorexit("Stencil not known to CFDSolver.");
	}
	if (m_stencil->getD() != dim) {
		natrium_errorexit("Stencil has wrong dimension.");
	}

	/// Build streaming data object
	if (SEDG == configuration->getAdvectionScheme()) {
		// start timer
		TimerOutput::Scope timer_section(Timing::getTimer(), "System assembly");

		// if this string is "": matrix is reassembled and not read from file
		string whereAreTheStoredMatrices;
		if (configuration->isRestartAtLastCheckpoint()) {
			whereAreTheStoredMatrices = configuration->getOutputDirectory();
		}
		// create SEDG MinLee by reading the system matrices from files or assembling
		// TODO estimate time for assembly
		try {
			m_advectionOperator = boost::make_shared<SEDGMinLee<dim> >(
					m_problemDescription->getMesh(),
					m_problemDescription->getBoundaries(),
					configuration->getSedgOrderOfFiniteElement(), m_stencil,
					whereAreTheStoredMatrices,
					(CENTRAL == configuration->getSedgFluxType()));
		} catch (AdvectionSolverException & e) {
			natrium_errorexit(e.what());
		}
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
	m_i = m_iterationStart;

/// Calculate relaxation parameter and build collision model
	double tau = 0.0;
	double gamma = -1.0;
	double G = -5;
	if (BGK_STANDARD == configuration->getCollisionScheme()) {
		tau = BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
		m_collisionModel = boost::make_shared<BGKStandard>(tau,
				m_configuration->getTimeStepSize(), m_stencil);
	} else if (BGK_STEADY_STATE == configuration->getCollisionScheme()) {
		gamma = configuration->getBGKSteadyStateGamma();
		tau = BGKSteadyState::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil, gamma);
		m_collisionModel = boost::make_shared<BGKSteadyState>(tau,
				m_configuration->getTimeStepSize(), m_stencil, gamma);
	} else if (BGK_STANDARD_TRANSFORMED
			== configuration->getCollisionScheme()) {
		tau = BGKStandardTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
		m_collisionModel = boost::make_shared<BGKStandardTransformed>(tau,
				m_configuration->getTimeStepSize(), m_stencil);
	} else if (BGK_MULTIPHASE == configuration->getCollisionScheme()) {
		tau = BGKStandardTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
		G = m_configuration->getBGKPseudopotentialG();
		boost::shared_ptr<BGKPseudopotential<dim> > coll_tmp =
				boost::make_shared<BGKPseudopotential<dim> >(tau,
						m_configuration->getTimeStepSize(), m_stencil, G);
		coll_tmp->setAdvectionOperator(m_advectionOperator);
		m_collisionModel = coll_tmp;
	} else if (BGK_INCOMPRESSIBLE == configuration->getCollisionScheme()) {
		tau = BGKIncompressible::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
		m_collisionModel = boost::make_shared<BGKIncompressible>(tau,
				m_configuration->getTimeStepSize(), m_stencil);
	} else if (MRT_STANDARD == configuration->getCollisionScheme()) {
		tau = BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
		m_collisionModel = boost::make_shared<MRTStandard>(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), m_stencil);
	} else if (KBC_STANDARD == configuration->getCollisionScheme()) {
		tau = KBCStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
		m_collisionModel = boost::make_shared<KBCStandard>(tau,
				m_configuration->getTimeStepSize(), m_stencil);
	}

// initialize macroscopic variables
#ifdef WITH_TRILINOS_MPI
	map<dealii::types::global_dof_index, dealii::Point<dim> > supportPoints;
	m_advectionOperator->mapDoFsToSupportPoints(supportPoints);
	m_density.reinit(m_advectionOperator->getLocallyOwnedDofs(),
	MPI_COMM_WORLD);
	m_tmpDensity.reinit(m_advectionOperator->getLocallyOwnedDofs(),
	MPI_COMM_WORLD);
	for (size_t i = 0; i < dim; i++) {
		distributed_vector vi(m_advectionOperator->getLocallyOwnedDofs(),
		MPI_COMM_WORLD);
		m_velocity.push_back(vi);
		m_tmpVelocity.push_back(vi);
	}
#else
	size_t numberOfDoFs = this->getNumberOfDoFs();
	map<dealii::types::global_dof_index, dealii::Point<dim> > supportPoints;
	m_advectionOperator->mapDoFsToSupportPoints(supportPoints);
	m_density.reinit(numberOfDoFs);
	m_tmpDensity.reinit(numberOfDoFs);
	for (size_t i = 0; i < dim; i++) {
		distributed_vector vi(numberOfDoFs);
		m_velocity.push_back(vi);
		m_tmpVelocity.push_back(vi);
	}
#endif
	applyInitialDensities(m_density, supportPoints);
	applyInitialVelocities(m_velocity, supportPoints);
	m_residuumDensity = 1.0;
	m_residuumVelocity = 1.0;

	//distribution functions
#ifdef WITH_TRILINOS_MPI
	m_f.reinit(m_stencil->getQ(), m_advectionOperator->getLocallyOwnedDofs(),
	MPI_COMM_WORLD);
#else
	m_f.reinit(m_stencil->getQ(), m_advectionOperator->getNumberOfDoFs());
#endif

	/// Build time integrator
	if (RUNGE_KUTTA_5STAGE == configuration->getTimeIntegrator()) {
		m_timeIntegrator = boost::make_shared<
				RungeKutta5LowStorage<distributed_sparse_block_matrix,
						distributed_block_vector> >(
				configuration->getTimeStepSize(), m_f.getFStream());
	} else if (THETA_METHOD == configuration->getTimeIntegrator()) {
		m_timeIntegrator = boost::make_shared<
				ThetaMethod<distributed_sparse_block_matrix,
						distributed_block_vector> >(
				configuration->getTimeStepSize(), m_f.getFStream(),
				configuration->getThetaMethodTheta());
	} else if (EXPONENTIAL == configuration->getTimeIntegrator()) {
		m_timeIntegrator = boost::make_shared<
				ExponentialTimeIntegrator<distributed_sparse_block_matrix,
						distributed_block_vector> >(
				configuration->getTimeStepSize(), m_stencil->getQ() - 1);
	} else if (OTHER == configuration->getTimeIntegrator()) {
		if (configuration->getDealIntegrator() < 7) {
			m_timeIntegrator = boost::make_shared<
					DealIIWrapper<distributed_sparse_block_matrix,
							distributed_block_vector> >(
					configuration->getTimeStepSize(),
					configuration->getDealIntegrator(),
					configuration->getDealLinearSolver());
		} else if (configuration->getDealIntegrator() < 12) {
			m_timeIntegrator =
					boost::make_shared<
							DealIIWrapper<distributed_sparse_block_matrix,
									distributed_block_vector> >(
							configuration->getTimeStepSize(),
							configuration->getDealIntegrator(),
							configuration->getDealLinearSolver(),
							configuration->getEmbeddedDealIntegratorCoarsenParameter(),
							configuration->getEmbeddedDealIntegratorRefinementParameter(),
							configuration->getEmbeddedDealIntegratorMinimumTimeStep(),
							configuration->getEmbeddedDealIntegratorMaximumTimeStep(),
							configuration->getEmbeddedDealIntegratorRefinementTolerance(),
							configuration->getEmbeddedDealIntegratorCoarsenTolerance());
		};
	}

// OUTPUT

	double maxU = getMaxVelocityNorm();
	double charU = problemDescription->getCharacteristicVelocity();
	if (charU == 0.0) {
		charU = maxU;
	}
	double dx = CFDSolverUtilities::getMinimumVertexDistance<dim>(
			*problemDescription->getMesh());
	LOG(WELCOME) << "------ NATriuM solver ------" << endl;
	LOG(WELCOME) << "------ commit " << Info::getGitSha() << " ------" << endl;
	LOG(WELCOME) << "------ " << currentDateTime() << " ------" << endl;
	LOG(WELCOME) << "------ " << Info::getUserName() << " on "
			<< Info::getHostName() << " ------" << endl;
	LOG(WELCOME) << "viscosity:                "
			<< problemDescription->getViscosity() << " m^2/s" << endl;
	LOG(WELCOME) << "char. length:             "
			<< problemDescription->getCharacteristicLength() << " m" << endl;
	LOG(WELCOME) << "max |u_0|:                "
			<< maxU * problemDescription->getCharacteristicLength() << " m/s"
			<< endl;
	LOG(WELCOME) << "Reynolds number:          "
			<< (charU * problemDescription->getCharacteristicLength())
					/ problemDescription->getViscosity() << endl;
	double Ma = charU / m_stencil->getSpeedOfSound();
	LOG(WELCOME) << "Mach number:              " << Ma << endl;
	LOG(WELCOME) << "Stencil scaling:          "
			<< configuration->getStencilScaling() << endl;
	LOG(WELCOME) << "Sound speed:              " << m_stencil->getSpeedOfSound()
			<< endl;
	//TODO propose optimal cfl based on time integrator
	const double optimal_cfl = 0.4;
	LOG(WELCOME) << "Recommended dt (CFL 0.4): "
			<< CFDSolverUtilities::calculateTimestep<dim>(
					*m_problemDescription->getMesh(),
					configuration->getSedgOrderOfFiniteElement(), *m_stencil,
					optimal_cfl) << " s" << endl;
	LOG(WELCOME) << "Actual dt:                "
			<< configuration->getTimeStepSize() << " s" << endl;
	LOG(WELCOME) << "CFL number:               "
			<< configuration->getTimeStepSize() / dx
					* m_stencil->getMaxParticleVelocityMagnitude()
					* configuration->getSedgOrderOfFiniteElement()
					* configuration->getSedgOrderOfFiniteElement() << endl;
	LOG(WELCOME) << "dx:                       " << dx << endl;
	LOG(WELCOME) << "Order of finite element:  "
			<< configuration->getSedgOrderOfFiniteElement() << endl;
	LOG(WELCOME) << "----------------------------" << endl;
	LOG(WELCOME) << "== COLLISION ==          " << endl;
	switch (configuration->getCollisionScheme()) {
	case BGK_STANDARD: {
		LOG(WELCOME) << "tau:                      " << tau << endl;
		break;
	}
	case BGK_STANDARD_TRANSFORMED: {
		LOG(WELCOME) << "tau:                      " << tau << endl;
		break;
	}
	case BGK_STEADY_STATE: {
		LOG(WELCOME) << "tau:                      " << tau << endl;
		LOG(WELCOME) << "steady state gamma:       " << gamma << endl;
		LOG(WELCOME) << "Effective Ma:             " << Ma / sqrt(gamma)
				<< endl;

		break;
	}
	case BGK_MULTIPHASE: {
		LOG(WELCOME) << "tau:                      " << tau << endl;
		LOG(WELCOME) << "G:                        " << G << endl;
		break;
	}
	case BGK_INCOMPRESSIBLE: {
		LOG(WELCOME) << "tau:                      " << tau << endl;
		break;
	}
	case MRT_STANDARD: {
		LOG(WELCOME) << "tau:						" << tau << endl;
		break;
	}
	case KBC_STANDARD: {
		LOG(WELCOME) << "tau:						" << tau << endl;
		break;
	}
	}
	LOG(WELCOME) << "----------------------------" << endl;

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
	DealIIExtensions::extract_dofs_with_support_on_boundary(
			*(m_advectionOperator->getDoFHandler()), dealii::ComponentMask(),
			m_isDoFAtBoundary, boundaryIndicators);
	size_t nofBoundaryNodes = 0;
	for (size_t i = 0; i < m_isDoFAtBoundary.size(); i++) {
		if (m_isDoFAtBoundary.at(i)) {
			nofBoundaryNodes += 1;
		}
	}
	LOG(DETAILED) << "Number of non-periodic boundary grid points: "
			<< nofBoundaryNodes << endl;
	LOG(DETAILED) << "Number of total grid points: " << getNumberOfDoFs()
			<< endl;
	const vector<dealii::types::global_dof_index>& dofs_per_proc =
			m_advectionOperator->getDoFHandler()->n_locally_owned_dofs_per_processor();
	for (size_t i = 0;
			i < dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); i++) {
		LOG(DETAILED) << "Process "
				<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
				<< " has " << dofs_per_proc.at(i) << " grid points." << endl;
	}

// Initialize distribution functions
	if (configuration->isRestartAtLastCheckpoint()) {
		loadDistributionFunctionsFromFiles(
				m_configuration->getOutputDirectory());
	} else {
		initializeDistributions();
	}
	// initialize nonlinear boundaries
	m_nonlinearBoundaryVector.reinit(m_f.getFStream());
	m_problemDescription->getBoundaries()->initializeNonlinearBoundaries(
			m_advectionOperator, m_stencil, &m_density, &m_velocity, &m_f,
			&m_nonlinearBoundaryVector);

// Create file for output table
	if ((not configuration->isSwitchOutputOff())
	/*and (configuration->getOutputTableInterval()
	 < configuration->getNumberOfTimeSteps())*/) {
		std::stringstream s;
		s << configuration->getOutputDirectory().c_str()
				<< "/results_table.txt";
		//create the SolverStats object which is responsible for the results table
		m_solverStats = boost::make_shared<SolverStats<dim> >(this, s.str());
	} else {
		m_solverStats = boost::make_shared<SolverStats<dim> >(this);
	}

}
/* Constructor */
template CFDSolver<2>::CFDSolver(
		boost::shared_ptr<SolverConfiguration> configuration,
		boost::shared_ptr<ProblemDescription<2> > problemDescription);
template CFDSolver<3>::CFDSolver(
		boost::shared_ptr<SolverConfiguration> configuration,
		boost::shared_ptr<ProblemDescription<3> > problemDescription);

template<size_t dim>
void CFDSolver<dim>::stream() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Streaming");

// no streaming in direction 0; begin with 1
	distributed_block_vector& f = m_f.getFStream();
	const distributed_sparse_block_matrix& systemMatrix =
			m_advectionOperator->getSystemMatrix();
	const distributed_block_vector& systemVector =
			m_advectionOperator->getSystemVector();

	try {
		m_time = m_timeIntegrator->step(f, systemMatrix, systemVector, m_time,
				m_timeIntegrator->getTimeStepSize());
	} catch (std::exception& e) {
		natrium_errorexit(e.what());
	}
	m_collisionModel->setTimeStep(m_timeIntegrator->getTimeStepSize());
}
template void CFDSolver<2>::stream();
template void CFDSolver<3>::stream();

template<size_t dim>
void CFDSolver<dim>::collide() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Collision");

	try {
		m_collisionModel->collideAll(m_f, m_density, m_velocity,
				m_advectionOperator->getLocallyOwnedDofs(), false);
	} catch (CollisionException& e) {
		natrium_errorexit(e.what());
	}
}
template void CFDSolver<2>::collide();
template void CFDSolver<3>::collide();

template<size_t dim>
void CFDSolver<dim>::reassemble() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Reassemble");

	try {
		m_advectionOperator->reassemble();
	} catch (AdvectionSolverException & e) {
		natrium_errorexit(e.what());
	}
}
template void CFDSolver<2>::reassemble();
template void CFDSolver<3>::reassemble();

template<size_t dim>
void CFDSolver<dim>::run() {
	m_i = m_iterationStart;
	while (true) {
		if (stopConditionMet()) {
			break;
		}
		output(m_i);
		m_i++;
		stream();
		collide();
	}
	output(m_i);

// Finalize
	Timing::getTimer().print_summary();
	LOG(BASIC) << "NATriuM run complete." << endl;
	LOG(BASIC) << "Summary: " << endl;
	LOG(BASIC) << Timing::getOutStream().str() << endl;
}
template void CFDSolver<2>::run();
template void CFDSolver<3>::run();

template<size_t dim>
bool CFDSolver<dim>::stopConditionMet() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(),
			"Check stop condition");

// Maximum number of iterations
	size_t N = m_configuration->getNumberOfTimeSteps();
	if (m_i >= N) {
		LOG(BASIC)
				<< "Stop condition: Maximum number of iterations reached in iteration "
				<< m_i << "." << endl;
		return true;
	}
// End time
	const double end_time = m_configuration->getSimulationEndTime();
	if (m_time >= end_time) {
		LOG(BASIC) << "Stop condition: Simulation end time t_max=" << end_time
				<< " reached in iteration " << m_i << "." << endl;
		return true;
	}
// Converged
	const size_t check_interval = 10;
	const double convergence_threshold =
			m_configuration->getConvergenceThreshold();
	if (m_i % 10 == 0) {
		m_solverStats->calulateResiduals(m_i);
		if ((m_residuumVelocity < convergence_threshold)
		/*and (m_residuumDensity < convergence_threshold)*/) {
			LOG(BASIC)
					<< "Stop condition: Simulation converged below threshold "
					<< convergence_threshold << " in iteration " << m_i << "."
					<< endl;
			LOG(BASIC) << "The actual variation was " << m_residuumVelocity
					<< " on velocity and " << m_residuumDensity
					<< " on density between iterations " << m_i - check_interval
					<< " and " << m_i << "." << endl;
			return true;
		}
	}
	return false;
}
template bool CFDSolver<2>::stopConditionMet();
template bool CFDSolver<3>::stopConditionMet();

template<size_t dim>
void CFDSolver<dim>::output(size_t iteration) {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Output");

// output: vector fields as .vtu files
	if (iteration == m_iterationStart) {
		m_tstart = time(0);
	}
	if (not m_configuration->isSwitchOutputOff()) {
		if (iteration % 100 == 0) {
			LOG(DETAILED) << "Iteration " << iteration << ",  t = " << m_time
					<< endl;
		}
		// output estimated runtime after iterations 1, 10, 100, 1000, ...
		/*if (iteration > m_iterationStart) {
		 if (int(log10(iteration - m_iterationStart))
		 == log10(iteration - m_iterationStart)) {
		 time_t estimated_end = m_tstart
		 + (m_configuration->getNumberOfTimeSteps()
		 - m_iterationStart)
		 / (iteration - m_iterationStart)
		 * (time(0) - m_tstart);
		 struct tm * ltm = localtime(&estimated_end);
		 LOG(BASIC) << "i = " << iteration << "; Estimated end: "
		 << string(asctime(ltm)) << endl;
		 }
		 }*/
		if (iteration % m_configuration->getOutputSolutionInterval() == 0) {
			// save local part of the solution
			std::stringstream str;
			str << m_configuration->getOutputDirectory().c_str() << "/t_"
					<< m_problemDescription->getMesh()->locally_owned_subdomain()
					<< "." << iteration << ".vtu";
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

			/// For Benchmarks: add analytic solution
			addAnalyticSolutionToOutput(data_out);

			// tell the data processor the locally owned cells
			dealii::Vector<float> subdomain(
					m_problemDescription->getMesh()->n_active_cells());
			for (unsigned int i = 0; i < subdomain.size(); ++i)
				subdomain(i) =
						m_problemDescription->getMesh()->locally_owned_subdomain();
			data_out.add_data_vector(subdomain, "subdomain");

			// Write vtu file
			data_out.build_patches(
					m_configuration->getSedgOrderOfFiniteElement() + 1);
			data_out.write_vtu(vtu_output);

			// Write pvtu file (which is a master file for all the single vtu files)
			if (is_MPI_rank_0()) {
				// generate .pvtu filename
				std::stringstream pvtu_filename;
				pvtu_filename << m_configuration->getOutputDirectory().c_str()
						<< "/t_"
						<< m_problemDescription->getMesh()->locally_owned_subdomain()
						<< "." << iteration << ".pvtu";
				std::ofstream pvtu_output(pvtu_filename.str().c_str());

				// generate all other filenames
				std::vector<std::string> filenames;
				for (unsigned int i = 0;
						i < dealii::Utilities::MPI::n_mpi_processes(
						MPI_COMM_WORLD); ++i) {
					std::stringstream vtu_filename_i;
					vtu_filename_i
							<< m_configuration->getOutputDirectory().c_str()
							<< "/t_" << i << "." << iteration << ".vtu";
					filenames.push_back(vtu_filename_i.str());
				}
				data_out.write_pvtu_record(pvtu_output, filenames);
			}
		}

		// output: table
		// calculate information + physical properties

		if (iteration % m_configuration->getOutputTableInterval() == 0) {
			m_solverStats->printNewLine();
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
// PRECONDITION: vectors already created with the right sizes

	LOG(BASIC) << "Initialize distribution functions: ";
	vector<double> feq(m_stencil->getQ());
	numeric_vector u(dim);

// save starting time
	double t0 = m_time;

// Initialize f with the equilibrium distribution functions
//for all degrees of freedom on current processor
	const dealii::IndexSet& locally_owned_dofs =
			m_advectionOperator->getLocallyOwnedDofs();
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;
		for (size_t j = 0; j < dim; j++) {
			u(j) = m_velocity.at(j)(i);
		}
		m_collisionModel->getEquilibriumDistributions(feq, u, m_density(i));
		for (size_t j = 0; j < m_stencil->getQ(); j++) {
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
		double residual = 1000000000;
		const bool inInitializationProcedure = true;
		distributed_vector oldDensities;
		while (residual > m_configuration->getIterativeInitializationResidual()) {
			if (loopCount
					> m_configuration->getIterativeInitializationNumberOfIterations()) {
				LOG(WARNING)
						<< "The iterative Initialization of equilibrium distribution functions could only reach residual "
						<< residual << " after " << loopCount
						<< " iterations (Aimed at residual "
						<< m_configuration->getIterativeInitializationResidual()
						<< "). If that is too bad, increase the number of iterations in the iterative initialization scheme. "
						<< "To avoid this Warning, soften the scheme (i.e. aim at a greater residual.)";
				break;
			}
			oldDensities = m_density;
			try {
				stream();
			} catch (std::exception& e) {
				natrium_errorexit(e.what());
			}
			// collide without recalculating velocities
			try {
				m_collisionModel->collideAll(m_f, m_density, m_velocity,
						m_advectionOperator->getLocallyOwnedDofs(),
						inInitializationProcedure);
			} catch (CollisionException& e) {
				natrium_errorexit(e.what());
			}
			oldDensities -= m_density;
			residual = oldDensities.norm_sqr();
			loopCount++;
		}
		LOG(DETAILED) << "Residual " << residual << " reached after "
				<< loopCount << " iterations." << endl;

		//for all degrees of freedom on current processor
		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;
			for (size_t j = 0; j < dim; j++) {
				u(j) = m_velocity.at(j)(i);
			}
			m_collisionModel->getEquilibriumDistributions(feq, u, m_density(i));
			for (size_t j = 0; j < m_stencil->getQ(); j++) {
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

	m_time = t0;

	LOG(BASIC) << "Initialize distribution functions: done." << endl;
}

template void CFDSolver<2>::initializeDistributions();
template void CFDSolver<3>::initializeDistributions();

template<size_t dim>
void CFDSolver<dim>::saveDistributionFunctionsToFiles(const string& directory) {
	for (size_t i = 0; i < m_stencil->getQ(); i++) {
		// filename
		std::stringstream filename;

		filename << directory << "/checkpoint_f_" << i
#ifdef WITH_TRILINOS_MPI
				<< "."
				<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
#endif
				<< ".dat";
		std::ofstream file(filename.str().c_str());
#ifndef WITH_TRILINOS
		m_f.at(i).block_write(file);
#else
		// TODO Write and read functions for Trilinos vectors. This here is really bad.
		numeric_vector tmp(m_f.at(i));
		tmp.block_write(file);
#endif
	}
}
template void CFDSolver<2>::saveDistributionFunctionsToFiles(
		const string& directory);
template void CFDSolver<3>::saveDistributionFunctionsToFiles(
		const string& directory);

template<size_t dim>
void CFDSolver<dim>::loadDistributionFunctionsFromFiles(
		const string& directory) {
// PRECONDITION: vectors already created with the right sizes
// read the distribution functions from file
	try {
		for (size_t i = 0; i < m_stencil->getQ(); i++) {
			// filename
			std::stringstream filename;
			filename << directory << "/checkpoint_f_" << i
#ifdef WITH_TRILINOS_MPI
					<< "."
					<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
#endif
					<< ".dat";
			std::ifstream file(filename.str().c_str());
#ifndef WITH_TRILINOS
			m_f.at(i).block_read(file);
#else
			// TODO Write and read functions for Trilinos vectors. This here is really bad.
			numeric_vector tmp(m_f.at(i));
			tmp.block_read(file);
			m_f.at(i) = tmp;

#endif
		}
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		natrium_errorexit(
				"An error occurred while reading the distribution functions from file: Please switch off the restart option to start the simulation from the beginning.");
	}
}
template void CFDSolver<2>::loadDistributionFunctionsFromFiles(
		const string& directory);
template void CFDSolver<3>::loadDistributionFunctionsFromFiles(
		const string& directory);

template<size_t dim>
void natrium::CFDSolver<dim>::applyInitialDensities(
		distributed_vector& initialDensities,
		const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
// get Function instance
	const boost::shared_ptr<dealii::Function<dim> >& f_rho =
			m_problemDescription->getInitialRhoFunction();
	const unsigned int dofs_per_cell =
			m_advectionOperator->getFe()->dofs_per_cell;
	vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			m_advectionOperator->getDoFHandler()->begin_active(), endc =
			m_advectionOperator->getDoFHandler()->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			cell->get_dof_indices(local_dof_indices);
			for (size_t i = 0; i < dofs_per_cell; i++) {
				assert(
						initialDensities.in_local_range(
								local_dof_indices.at(i)));
				assert(
						supportPoints.find(local_dof_indices.at(i))
								!= supportPoints.end());
				initialDensities(local_dof_indices.at(i)) = f_rho->value(
						supportPoints.at(local_dof_indices.at(i)));
			}
		} /* if is locally owned */
	} /* for all cells */
}
template
void natrium::CFDSolver<2>::applyInitialDensities(
		distributed_vector& initialDensities,
		const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints) const;
template
void natrium::CFDSolver<3>::applyInitialDensities(
		distributed_vector& initialDensities,
		const map<dealii::types::global_dof_index, dealii::Point<3> >& supportPoints) const;

template<size_t dim>
void natrium::CFDSolver<dim>::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
// get Function instance
	const boost::shared_ptr<dealii::Function<dim> >& f_u =
			m_problemDescription->getInitialUFunction();
	const unsigned int dofs_per_cell =
			m_advectionOperator->getFe()->dofs_per_cell;
	vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			m_advectionOperator->getDoFHandler()->begin_active(), endc =
			m_advectionOperator->getDoFHandler()->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			cell->get_dof_indices(local_dof_indices);
			for (size_t i = 0; i < dofs_per_cell; i++) {
				assert(
						initialVelocities.at(0).in_local_range(
								local_dof_indices.at(i)));
				assert(
						initialVelocities.at(1).in_local_range(
								local_dof_indices.at(i)));
				assert(
						supportPoints.find(local_dof_indices.at(i))
								!= supportPoints.end());
				for (size_t component = 0; component < dim; component++) {
					initialVelocities.at(component)(local_dof_indices.at(i)) =
							f_u->value(
									supportPoints.at(local_dof_indices.at(i)),
									component);
				}
			}
		} /* if is locally owned */
	} /* for all cells */
}
template
void natrium::CFDSolver<2>::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints) const;
template
void natrium::CFDSolver<3>::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const map<dealii::types::global_dof_index, dealii::Point<3> >& supportPoints) const;

template<size_t dim>
double CFDSolver<dim>::getTau() const {
	if ((BGK_STANDARD == m_configuration->getCollisionScheme())
			or (BGK_STANDARD_TRANSFORMED
					== m_configuration->getCollisionScheme())) {
		return BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(),
				m_configuration->getTimeStepSize(), *m_stencil);
	}
	LOG(WARNING)
			<< "getTau() is called, but you don't have a BGK Standard model."
			<< endl;
	return 0;
}
template double CFDSolver<2>::getTau() const;
template double CFDSolver<3>::getTau() const;

} /* namespace natrium */

