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
#include "Checkpoint.h"

#include "../stencils/D2Q9.h"
#include "../stencils/D3Q19.h"
#include "../stencils/D3Q15.h"
#include "../stencils/D3Q27.h"
#include "../stencils/Stencil.h"

#include "../advection/SEDGMinLee.h"
#include "../advection/SemiLagrangian.h"

#include "../collision/BGKStandard.h"
#include "../collision/BGKStandardTransformed.h"
#include "../collision/BGKSteadyState.h"
#include "../collision/BGKIncompressible.h"
#include "../collision/MRTStandard.h"
#include "../collision/KBCStandard.h"
#include "../collision/KBCCentral.h"
#include "../collision/BGKMultistep.h"

#include "../smoothing/ExponentialFilter.h"
#include "../smoothing/NewFilter.h"

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

	boost::filesystem::path out_dir(configuration->getOutputDirectory());
	boost::filesystem::path log_file = out_dir / "natrium.log";
	boost::filesystem::path checkpoint_dir = out_dir / "checkpoint";

	// CONFIGURE LOGGER
	if (configuration->isSwitchOutputOff()) {
		LOGGER().setConsoleLevel(SILENT);
		LOGGER().setFileLevel(SILENT);
	} else {
		LOGGER().setConsoleLevel(
				LogLevel(configuration->getCommandLineVerbosity()));
		LOGGER().setFileLevel(ALL);
		LOGGER().setLogFile(log_file.string());
	}

	// welcome message
	LOG(WELCOME) << "------ NATriuM solver ------" << endl;
	LOG(WELCOME) << "------ commit " << Info::getGitSha() << " ------" << endl;
	LOG(WELCOME) << "------ " << currentDateTime() << " ------" << endl;
	LOG(WELCOME) << "------ " << Info::getUserName() << " on "
			<< Info::getHostName() << " ------" << endl;

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

	/// Build stencil
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

	// Restart from checkpoint?
	boost::shared_ptr<Checkpoint<dim> > checkpoint;
	CheckpointStatus checkpoint_status;
	size_t restart_i = m_configuration->getRestartAtIteration();
	if (0 != restart_i) {
		checkpoint = boost::make_shared<Checkpoint<dim> >(restart_i,
				checkpoint_dir);
		if (not checkpoint->exists()) {
			std::stringstream msg;
			msg << "You want to restart from iteration " << restart_i
					<< ", but I could not find the required checkpoint files "
					<< checkpoint->getStatusFile().string() << " and "
					<< checkpoint->getDataFile().string() << ".";
			natrium_errorexit(msg.str().c_str());
		} else {
			LOG(BASIC) << "Restart at iteration " << restart_i << endl;
		}
	}

	/// Build streaming data object
	LOG(WELCOME) << "Create streaming object ..." << endl;
	if (SEDG == configuration->getAdvectionScheme()) {
		// start timer

		// create SEDG MinLee and assemble
		try {
			m_advectionOperator = boost::make_shared<SEDGMinLee<dim> >(
					m_problemDescription->getMesh(),
					m_problemDescription->getBoundaries(),
					configuration->getSedgOrderOfFiniteElement(), m_stencil,
					(CENTRAL == configuration->getSedgFluxType()));
		} catch (AdvectionSolverException & e) {
			natrium_errorexit(e.what());
		}
	} else if (SEMI_LAGRANGIAN == configuration->getAdvectionScheme()) {
		try {
			m_advectionOperator = boost::make_shared<SemiLagrangian<dim> >(
					m_problemDescription->getMesh(),
					m_problemDescription->getBoundaries(),
					configuration->getSedgOrderOfFiniteElement(), m_stencil,
					0.0);
		} catch (AdvectionSolverException & e) {
			natrium_errorexit(e.what());
		}
	}
	// Refine mesh, Build DoF system (by reading from restart file)
	if (checkpoint) {
		checkpoint->load(m_f, *m_problemDescription, *m_advectionOperator,
				checkpoint_status);
		m_iterationStart = checkpoint_status.iterationNumber;
		m_time = checkpoint_status.time;
	} else {
		m_problemDescription->refineAndTransform();
		m_advectionOperator->setupDoFs();
		m_f.reinit(m_stencil->getQ(),
				m_advectionOperator->getLocallyOwnedDofs(),
				MPI_COMM_WORLD);
		m_iterationStart = 0;
		m_time = 0;
	}
	m_i = m_iterationStart;

	// set time step size
	double delta_t = CFDSolverUtilities::calculateTimestep<dim>(
			*(m_problemDescription->getMesh()),
			m_configuration->getSedgOrderOfFiniteElement(), *m_stencil,
			configuration->getCFL());

	// assemble advection operator
	// upon setDeltaT, the semi-Lagrangian updates the sparsity pattern
	{
		TimerOutput::Scope timer_section(Timing::getTimer(), "Assembly");
		m_advectionOperator->setDeltaT(delta_t);
		m_advectionOperator->reassemble();
	}
	LOG(WELCOME) << "... done" << endl;

/// Calculate relaxation parameter and build collision model
	double tau = 0.0;
	double gamma = -1.0;
	double G = -5;

	if (BGK_STANDARD == configuration->getCollisionScheme()) {
		tau = BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		m_collisionModel = boost::make_shared<BGKStandard>(tau, delta_t,
				m_stencil);
	} else if (BGK_STEADY_STATE == configuration->getCollisionScheme()) {
		gamma = configuration->getBGKSteadyStateGamma();
		tau = BGKSteadyState::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil,
				gamma);
		m_collisionModel = boost::make_shared<BGKSteadyState>(tau, delta_t,
				m_stencil, gamma);
	} else if (BGK_STANDARD_TRANSFORMED
			== configuration->getCollisionScheme()) {
		tau = BGKStandardTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		m_collisionModel = boost::make_shared<BGKStandardTransformed>(tau,
				delta_t, m_stencil);
	} else if (BGK_MULTIPHASE == configuration->getCollisionScheme()) {
		tau = BGKStandardTransformed::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		PseudopotentialParameters pp_para(
				m_configuration->getPseudopotentialType(),
				m_configuration->getBGKPseudopotentialG(),
				m_configuration->getBGKPseudopotentialT());
		G = m_configuration->getBGKPseudopotentialG();
		boost::shared_ptr<BGKPseudopotential<dim> > coll_tmp =
				boost::make_shared<BGKPseudopotential<dim> >(tau, delta_t,
						m_stencil, pp_para);
		coll_tmp->setAdvectionOperator(m_advectionOperator);
		m_collisionModel = coll_tmp;
	} else if (BGK_INCOMPRESSIBLE == configuration->getCollisionScheme()) {
		tau = BGKIncompressible::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		m_collisionModel = boost::make_shared<BGKIncompressible>(tau, delta_t,
				m_stencil);
	} else if (MRT_STANDARD == configuration->getCollisionScheme()) {
		tau = BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		m_collisionModel = boost::make_shared<MRTStandard>(
				m_problemDescription->getViscosity(), delta_t, m_stencil);
	} else if (KBC_STANDARD == configuration->getCollisionScheme()) {
		tau = KBCStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		m_collisionModel = boost::make_shared<KBCStandard>(tau, delta_t,
				m_stencil);
	} else if (KBC_CENTRAL == configuration->getCollisionScheme()) {
		tau = KBCCentral::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		m_collisionModel = boost::make_shared<KBCCentral>(tau, delta_t,
				m_stencil);
	} else if (BGK_MULTI_AM4 == configuration->getCollisionScheme()) {
		tau = BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		boost::shared_ptr<BGKMultistep> bgk_tmp = boost::make_shared<
				BGKMultistep>(tau, delta_t, m_stencil, 0); //0 for Adams Moulton
		m_multistepData = bgk_tmp;
		m_collisionModel = bgk_tmp;
		m_collisionModel->setViscosity(m_problemDescription->getViscosity());
	} else if (BGK_MULTI_BDF2 == configuration->getCollisionScheme()) {
		tau = BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
		boost::shared_ptr<BGKMultistep> bgk_tmp = boost::make_shared<
				BGKMultistep>(tau, delta_t, m_stencil, 1); // 1 for BDF2

		m_multistepData = bgk_tmp;
		m_collisionModel = bgk_tmp;
		m_collisionModel->setViscosity(m_problemDescription->getViscosity());
		// TODO call setViscosity only once (for general collision model)
		// TODO remove relaxation-parameter from collision model and calculate in each time step from dt and nu
}

	// apply external force
	if (m_problemDescription->hasExternalForce()) {
		m_collisionModel->setForceType(m_configuration->getForcingScheme());
		m_collisionModel->setExternalForce(
				*(m_problemDescription->getExternalForce()));
	}

// initialize macroscopic variables
	m_advectionOperator->mapDoFsToSupportPoints(m_supportPoints);
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
	if (not checkpoint) {
		applyInitialDensities(m_density, m_supportPoints);
		applyInitialVelocities(m_velocity, m_supportPoints);
	} else {
		calculateDensitiesAndVelocities();
	}
	m_residuumDensity = 1.0;
	m_residuumVelocity = 1.0;

	/// Build time integrator
	if (RUNGE_KUTTA_5STAGE == configuration->getTimeIntegrator()) {
		m_timeIntegrator = boost::make_shared<
				RungeKutta5LowStorage<distributed_sparse_block_matrix,
						distributed_block_vector> >(delta_t, m_f.getFStream());
	} else if (THETA_METHOD == configuration->getTimeIntegrator()) {
		m_timeIntegrator = boost::make_shared<
				ThetaMethod<distributed_sparse_block_matrix,
						distributed_block_vector> >(delta_t, m_f.getFStream(),
				configuration->getThetaMethodTheta());
	} else if (EXPONENTIAL == configuration->getTimeIntegrator()) {
		m_timeIntegrator = boost::make_shared<
				ExponentialTimeIntegrator<distributed_sparse_block_matrix,
						distributed_block_vector> >(delta_t,
				m_stencil->getQ() - 1);
	} else if (OTHER == configuration->getTimeIntegrator()) {
		if (configuration->getDealIntegrator() < 7) {
			m_timeIntegrator = boost::make_shared<
					DealIIWrapper<distributed_sparse_block_matrix,
							distributed_block_vector> >(delta_t,
					configuration->getDealIntegrator(),
					configuration->getDealLinearSolver());
		} else if (configuration->getDealIntegrator() < 12) {
			double min_delta_t = CFDSolverUtilities::calculateTimestep<dim>(
					*(m_problemDescription->getMesh()),
					m_configuration->getSedgOrderOfFiniteElement(), *m_stencil,
					configuration->getEmbeddedDealIntegratorMinimumCFL());
			double max_delta_t = CFDSolverUtilities::calculateTimestep<dim>(
					*(m_problemDescription->getMesh()),
					m_configuration->getSedgOrderOfFiniteElement(), *m_stencil,
					configuration->getEmbeddedDealIntegratorMaximumCFL());
			m_timeIntegrator =
					boost::make_shared<
							DealIIWrapper<distributed_sparse_block_matrix,
									distributed_block_vector> >(delta_t,
							configuration->getDealIntegrator(),
							configuration->getDealLinearSolver(),
							configuration->getEmbeddedDealIntegratorCoarsenParameter(),
							configuration->getEmbeddedDealIntegratorRefinementParameter(),
							min_delta_t, max_delta_t,
							configuration->getEmbeddedDealIntegratorRefinementTolerance(),
							configuration->getEmbeddedDealIntegratorCoarsenTolerance());
		};
	}
	if (m_problemDescription->getBoundaries()->hasNonlinearBoundaries()) {
		m_timeIntegrator->setBoundaryCollection(
				m_problemDescription->getBoundaries());
	}

// build filter
	if (m_configuration->isFiltering() == true) {
		if (m_configuration->getFilteringScheme() == EXPONENTIAL_FILTER) {
			m_filter = boost::make_shared<ExponentialFilter<dim> >(
					m_configuration->getExponentialFilterAlpha(),
					m_configuration->getExponentialFilterS(),
					m_configuration->getExponentialFilterNc(),
					m_configuration->isFilterDegreeByComponentSums(),
					*m_advectionOperator->getQuadrature(),
					*m_advectionOperator->getFe());
		} else if (m_configuration->getFilteringScheme() == NEW_FILTER) {
			m_filter = boost::make_shared<NewFilter<dim> >(
					m_configuration->getExponentialFilterAlpha(),
					m_configuration->getExponentialFilterS(),
					*m_advectionOperator->getQuadrature(),
					*m_advectionOperator->getFe());
		}
	}

// OUTPUT

	double maxU = getMaxVelocityNorm();
	double charU = problemDescription->getCharacteristicVelocity();
	if (charU == 0.0) {
		charU = maxU;
	}
	double dx = CFDSolverUtilities::getMinimumVertexDistance<dim>(
			*problemDescription->getMesh());
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
	LOG(WELCOME) << "----------------------------" << endl;
	LOG(WELCOME) << "== ADVECTION ==          " << endl;
	if (SEMI_LAGRANGIAN == configuration->getAdvectionScheme()) {
		LOG(WELCOME) << "Semi-Lagrangian advection" << endl;
		LOG(WELCOME) << "Semi-Lagrangian advection" << endl;
	} else if (SEDG == configuration->getAdvectionScheme()) {
		LOG(WELCOME) << "Spectral-element discontinuous Galerkin" << endl;
		const double optimal_cfl = 0.4;
		LOG(WELCOME) << "Recommended dt (CFL 0.4): "
				<< CFDSolverUtilities::calculateTimestep<dim>(
						*m_problemDescription->getMesh(),
						configuration->getSedgOrderOfFiniteElement(),
						*m_stencil, optimal_cfl) << " s" << endl;
	}
	LOG(WELCOME) << "Actual dt:                " << delta_t << " s" << endl;
	LOG(WELCOME) << "CFL number:               " << configuration->getCFL()
			<< endl;
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
	case KBC_CENTRAL: {
		LOG(WELCOME) << "tau:						" << tau << endl;
		break;
	}
	case BGK_MULTI_AM4: {
		LOG(WELCOME) << "tau:						" << tau << endl;
		break;
	}
	case BGK_MULTI_BDF2: {
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
		LOG(DETAILED) << "Process " << i << " has " << dofs_per_proc.at(i)
				<< " grid points." << endl;
	}

// Initialize distribution functions
	if (not checkpoint) {
		initializeDistributions();
	}

// initialize nonlinear boundaries
	m_boundaryVector.reinit(m_advectionOperator->getSystemVector());
	m_boundaryVector = m_advectionOperator->getSystemVector();
	m_problemDescription->getBoundaries()->initializeNonlinearBoundaries(
			m_advectionOperator, m_stencil, &m_density, &m_velocity, &m_f,
			&m_boundaryVector);

// Create file for output table
	if ((not configuration->isSwitchOutputOff())
	/*and (configuration->getOutputTableInterval()
	 < configuration->getNumberOfTimeSteps())*/) {
		std::stringstream s;
		s << configuration->getOutputDirectory().c_str()
				<< "/results_table.txt";
		//create the SolverStats object which is responsible for the results table
		m_solverStats = boost::make_shared<SolverStats<dim> >(this, s.str());
		//create the TurbulenceStats object which is responsible for the turbulence table
		if (configuration->isOutputTurbulenceStatistics()) {
			std::stringstream s2;
			s2 << configuration->getOutputDirectory().c_str()
					<< "/turbulence_table.txt";
			m_turbulenceStats = boost::make_shared<TurbulenceStats<dim> >(this,
					m_configuration->getWallNormalDirection(),
					m_configuration->getWallNormalCoordinates(), s2.str());
		}
		if (configuration->isOutputGlobalTurbulenceStatistics()){
			appendDataProcessor(boost::make_shared<GlobalTurbulenceStats<dim> > (*this));
		}
	} else {
		m_solverStats = boost::make_shared<SolverStats<dim> >(this);
		//m_turbulenceStats = boost::make_shared<TurbulenceStats<dim> >(this);
	}

	// print out memory requirements of single components
	LOG(BASIC) << endl << " ------- Memory Requirements (MPI rank 0) -------- " << endl;
	LOG(BASIC) << " |  Sparse matrix        |  "
			<< m_advectionOperator->getSystemMatrix().memory_consumption()
			<< " (#nonzero elem: "
			<< m_advectionOperator->getSystemMatrix().n_nonzero_elements()
			<< ")" << endl;
	LOG(BASIC) << " |  Mesh                 |  "
			<< m_problemDescription->getMesh()->memory_consumption()  << endl;
	LOG(BASIC) << " |  Distributions        |  " << m_f.memory_consumption()
			<< endl;
	LOG(BASIC) << " |  Tmp distributions    |  " << m_f.memory_consumption()
			<< endl;
	LOG(BASIC) << " |  Velocities           |  "
			<< dim * m_velocity.at(0).memory_consumption() << endl;
	LOG(BASIC) << " |  Densities            |  " << m_density.memory_consumption()
			<< endl;
	LOG(BASIC) << " ------------------------------------------------- " << endl << endl;

}
/* Constructor */

template<size_t dim>
void CFDSolver<dim>::stream() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Streaming");

// no streaming in direction 0; begin with 1
	distributed_block_vector& f = m_f.getFStream();
	const distributed_sparse_block_matrix& systemMatrix =
			m_advectionOperator->getSystemMatrix();

//const distributed_block_vector& systemVector =
//		m_advectionOperator->getSystemVector();
//TODO has to be replaced by boundaryVector
	m_boundaryVector = m_advectionOperator->getSystemVector();

//try {
	if (SEMI_LAGRANGIAN == m_configuration->getAdvectionScheme()) {
		distributed_block_vector f_tmp(f.n_blocks());
		reinitVector(f_tmp, f);
		f_tmp = f;
		//f_tmp.reinit(f);
		systemMatrix.vmult(f, f_tmp);
		//f += m_boundaryVector;

		if ((BGK_MULTI_AM4 == m_configuration->getCollisionScheme()
				|| (BGK_MULTI_BDF2 == m_configuration->getCollisionScheme()))
				&& (m_i - m_iterationStart) > 1) {
			distributed_block_vector& formerF =
					m_multistepData->getFormerF().getFStream();
			f_tmp = formerF;
			assert(m_multistepData != NULL);
			systemMatrix.vmult(formerF, f_tmp);

			distributed_block_vector& formerFEq =
					m_multistepData->getFormerFEq().getFStream();
			f_tmp = formerFEq;
			systemMatrix.vmult(formerFEq, f_tmp);
		}

		m_time += getTimeStepSize();

	} else {
		double new_dt = m_timeIntegrator->step(f, systemMatrix,
				m_boundaryVector, 0.0, m_timeIntegrator->getTimeStepSize());

		// For multistep methods, the former PDF for f and feq have also to be streamed
		if ((BGK_MULTI_AM4 == m_configuration->getCollisionScheme()
				|| (BGK_MULTI_BDF2 == m_configuration->getCollisionScheme()))
				&& (m_i - m_iterationStart) > 1) {

			distributed_block_vector& formerF =
					m_multistepData->getFormerF().getFStream();
			distributed_block_vector& formerFEq =
					m_multistepData->getFormerFEq().getFStream();

			m_timeIntegrator->step(formerF, systemMatrix, m_boundaryVector, 0.0,
					m_timeIntegrator->getTimeStepSize());

			m_timeIntegrator->step(formerFEq, systemMatrix, m_boundaryVector,
					0.0, m_timeIntegrator->getTimeStepSize());

		}

		m_timeIntegrator->setTimeStepSize(new_dt);
		m_time += new_dt;
		m_collisionModel->setTimeStep(m_timeIntegrator->getTimeStepSize());
	}
//} catch (std::exception& e) {
//	natrium_errorexit(e.what());
//}
}

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

template<size_t dim>
void CFDSolver<dim>::reassemble() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Reassembly");

	try {
		m_advectionOperator->reassemble();
	} catch (AdvectionSolverException & e) {
		natrium_errorexit(e.what());
	}
}

template<size_t dim>
void CFDSolver<dim>::filter() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Filter");

	if (m_configuration->isFiltering()) {
		if (m_i % m_configuration->getFilterInterval() == 0) {
			for (size_t i = 0; i < m_stencil->getQ(); i++) {
				m_filter->applyFilter(*m_advectionOperator->getDoFHandler(),
						m_f.at(i));
			}
		}
	}
}

template<size_t dim>
void CFDSolver<dim>::run() {
	m_i = m_iterationStart;
	collide();
	while (true) {
		if (stopConditionMet()) {
			break;
		}
		output(m_i);
		m_i++;
		stream();
		filter();
		collide();
		for (size_t i = 0; i < m_dataProcessors.size(); i++) {
			m_dataProcessors.at(i)->apply();
		}
	}
	output(m_i, true);

// Finalize
	Timing::getTimer().print_summary();
	LOG(BASIC) << "NATriuM run complete." << endl;
	LOG(BASIC) << "Summary: " << endl;
	LOG(BASIC) << Timing::getOutStream().str() << endl;
}

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
	const size_t check_interval = 100;
	const double convergence_threshold =
			m_configuration->getConvergenceThreshold();
	if (m_i % check_interval == 0) {
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

template<size_t dim>
void CFDSolver<dim>::output(size_t iteration, bool is_final) {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Output");

// sync MPI processes
	MPI_sync();

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
		if (m_configuration->isOutputTurbulenceStatistics())
			m_turbulenceStats->addToReynoldsStatistics(m_velocity);
		if ((iteration % m_configuration->getOutputSolutionInterval() == 0) or is_final) {
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
			/// For turbulent flows: add turbulent statistics
			if (m_configuration->isOutputTurbulenceStatistics()) {
				m_turbulenceStats->addReynoldsStatisticsToOutput(data_out);
			}

			// tell the data processor the locally owned cells
			dealii::Vector<float> subdomain(
					m_problemDescription->getMesh()->n_active_cells());
			for (unsigned int i = 0; i < subdomain.size(); ++i)
				subdomain(i) =
						m_problemDescription->getMesh()->locally_owned_subdomain();
			data_out.add_data_vector(subdomain, "subdomain");

			// Write vtu file
			data_out.build_patches(
					m_configuration->getSedgOrderOfFiniteElement() * 2);
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
					//<< m_configuration->getOutputDirectory().c_str() << "/"
					<< "t_" << i << "." << iteration << ".vtu";
					filenames.push_back(vtu_filename_i.str());
				}
				data_out.write_pvtu_record(pvtu_output, filenames);
			}
		}

		// output: table
		// calculate information + physical properties
		if (iteration % m_configuration->getOutputTableInterval() == 0) {
			m_solverStats->printNewLine();
			if (m_configuration->isOutputTurbulenceStatistics()) {
				assert(m_turbulenceStats);
				m_turbulenceStats->printNewLine();
			}
		}

		// output: checkpoint
		if ((iteration % m_configuration->getOutputCheckpointInterval() == 0) or is_final) {
			boost::filesystem::path checkpoint_dir(
					m_configuration->getOutputDirectory());
			checkpoint_dir /= "checkpoint";
			Checkpoint<dim> checkpoint(m_i, checkpoint_dir);
			CheckpointStatus checkpoint_status;
			checkpoint_status.iterationNumber = m_i;
			checkpoint_status.stencilScaling = m_stencil->getScaling();
			checkpoint_status.time = m_time;
			checkpoint_status.feOrder =
					m_configuration->getSedgOrderOfFiniteElement();
			checkpoint.write(*m_problemDescription->getMesh(), m_f,
					*m_advectionOperator->getDoFHandler(), checkpoint_status);
		} /*if checkpoint interval*/
	} /*if not output off*/
}

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

template<size_t dim>
double CFDSolver<dim>::getTau() const {
	double delta_t = m_timeIntegrator->getTimeStepSize();
	if ((BGK_STANDARD == m_configuration->getCollisionScheme())
			or (BGK_STANDARD_TRANSFORMED
					== m_configuration->getCollisionScheme())) {
		return BGKStandard::calculateRelaxationParameter(
				m_problemDescription->getViscosity(), delta_t, *m_stencil);
	}
	LOG(WARNING)
			<< "getTau() is called, but you don't have a BGK Standard model."
			<< endl;
	return 0;
}

template<size_t dim>
void CFDSolver<dim>::addToVelocity(
		boost::shared_ptr<dealii::Function<dim> > function) {
// get Function instance
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
						m_velocity.at(0).in_local_range(
								local_dof_indices.at(i)));
				assert(
						m_velocity.at(1).in_local_range(
								local_dof_indices.at(i)));
				assert(
						m_supportPoints.find(local_dof_indices.at(i))
								!= m_supportPoints.end());
				for (size_t component = 0; component < dim; component++) {
					m_velocity.at(component)(local_dof_indices.at(i)) =
							m_velocity.at(component)(local_dof_indices.at(i))
									+ function->value(
											m_supportPoints.at(
													local_dof_indices.at(i)),
											component);
				}
			}
		} /* if is locally owned */
	} /* for all cells */
	initializeDistributions();
}

template<size_t dim>
void CFDSolver<dim>::scaleVelocity(double scaling_factor) {
	m_velocity.at(0) *= scaling_factor;
	m_velocity.at(1) *= scaling_factor;
	if (dim == 3) {
		m_velocity.at(2) *= scaling_factor;
	}
	initializeDistributions();
}

template<size_t dim>
void CFDSolver<dim>::calculateDensitiesAndVelocities() {
// inefficient version, only for initialization

	size_t Q = m_stencil->getQ();
	m_density = 0;
	for (size_t i = 0; i < dim; i++) {
		m_velocity.at(i) = 0;
	}

	for (size_t i = 0; i < Q; i++) {
		m_density.add(m_f.at(i));
		for (size_t j = 0; j < dim; j++) {
			m_velocity.at(j).add(m_stencil->getDirection(i)(j), m_f.at(i));
		}
	}

//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(
			m_advectionOperator->getLocallyOwnedDofs().begin());
	dealii::IndexSet::ElementIterator end(
			m_advectionOperator->getLocallyOwnedDofs().end());
	for (; it != end; it++) {
		size_t i = *it;
		double one_by_rho_i = 1.0 / m_density(i);
		for (size_t j = 0; j < dim; j++) {
			m_velocity[j](i) = m_velocity[j](i) * one_by_rho_i;
		}
	}
}

template<size_t dim>
void CFDSolver<dim>::convertDeprecatedCheckpoint() {
	// load
	boost::filesystem::path checkpoint_dir(
			m_configuration->getOutputDirectory());
	checkpoint_dir /= "checkpoint";
	CheckpointStatus checkpoint_status;
	Checkpoint<dim>::loadFromDeprecatedCheckpointVersion(m_f,
			*m_advectionOperator, checkpoint_dir.string(), checkpoint_status);

	// save
	Checkpoint<dim> checkpoint(checkpoint_status.iterationNumber,
			checkpoint_dir);
	checkpoint.write(*m_problemDescription->getMesh(), m_f,
			*m_advectionOperator->getDoFHandler(), checkpoint_status);
}

template class CFDSolver<2> ;
template class CFDSolver<3> ;

} /* namespace natrium */

