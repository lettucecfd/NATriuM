/**
 * @file SolverConfiguration.h
 * @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SOLVERCONFIGURATION_H_
#define SOLVERCONFIGURATION_H_

#include "deal.II/base/parameter_handler.h"

#include "../problemdescription/ProblemDescription.h"
#include "../boltzmannmodels/BoltzmannModel.h"
#include "../utilities/BasicNames.h"
#include "../utilities/Logging.h"

namespace natrium {

/**
 * @short Implemented streaming data types
 */
enum AdvectionOperatorType {
	Advection_SEDGMinLee
};

/**
 * @short Implemented collision models
 */
enum CollisionType {
	Collision_BGKTransformed // Collision for the transformed distribution function as defined in MinLee2011
};

// StencilType defined in BoltzmannModel.h

/**
 * @short Implemented time integrators
 */
enum TimeIntegratorType {
	Integrator_RungeKutta5LowStorage
};

/**
 * @short the numerical flux used to calculate the advection operator
 */
enum FluxType {
	Flux_LaxFriedrichs, Flux_Central
};

/**
 * Output flags
 */
enum OutputFlags {
	out_noOutput = 0,
	out_CommandLineError = 1,
	out_CommandLineBasic = 2,
	out_CommandLineFull = 4,
	out_LogFile = 8,
	out_VectorFields = 16,
	out_Checkpoints = 32
};

/**
 * @short the initialization procedure for the distribution functions
 */
enum DistributionInitType {
	Equilibrium, // Distribute with equilibrium functions
	Iterative // Distribute with iterative procedure; enforces consistent initial conditions
};

/**
 * @short Exception class for CFDSolver
 */
class ConfigurationException: public std::exception {
private:
	std::string message;
public:
	ConfigurationException(const char *msg) :
			message(msg) {
	}
	~ConfigurationException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/** @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 *  @tparam dim The dimension of the flow (2 or 3).
 */
class SolverConfiguration: public dealii::ParameterHandler {
private:

	/// Streaming data type (e.g. MinLee2011)
	AdvectionOperatorType m_advectionOperatorType;

	/// Collision type (e.g. BGKTransformed)
	CollisionType m_collisionType;

	/// Stencil type (e.g. D2Q9)
	StencilType m_stencilType;

	/// Time Integrator type (e.g. RK5LowStorage)
	TimeIntegratorType m_timeIntegratorType;

	/// Numerical Flux
	FluxType m_fluxType;

	/// Time step size
	double m_timeStep;

	/// Number of time steps
	size_t m_numberOfTimeSteps;

	/// Order of finite element
	size_t m_orderOfFiniteElement;

	/// scaling of the difference stencil for the discrete particle velocity
	double m_dQScaling;

	/// Output directory
	std::string m_outputDirectory;

	/// the output flags
	int m_outputFlags;

	/// output frequency
	size_t m_outputVectorFieldsEvery;
	size_t m_outputCheckpointEvery;

	/// restart option
	bool m_restart;

	/// initialization procedure
	DistributionInitType m_distributionInitType;

	/// max initialization iterations
	size_t m_maxDistributionInitIterations;

	/// stop condition for the density residual of the initialization procedure
	double m_stopDistributionInitResidual;

public:

	/// constructor
	SolverConfiguration() {
		// Declare structure of the parameter file
		enter_subsection("General");
		{
			declare_entry("Time step size", "0.2",
					dealii::Patterns::Double(1e-10),
					"Size of the (initial) time step.");
			declare_entry("Switch off output", "false",
					dealii::Patterns::Bool(), "Switch output off, completely.");
			declare_entry("Stencil", "D2Q9",
					dealii::Patterns::Selection("D2Q9"),
					"The discrete velocity stencil. The number behind D denotes the dimension (2 or 3). The number behind Q denotes the number of particle directions in the discrete velocity model.");
			declare_entry("Stencil scaling", "1.0",
					dealii::Patterns::Double(1e-10),
					"The scaling of the discrete velocities. Whereas in the standard LBM the magnitude of the particle velocities is set to 1.0 due to the uniform mesh grid, the SEDG LBM features scaled particle velocities. As the scaling factor is proportional to the speed of sound, it strongly impacts the relaxation time.");
		}
		leave_subsection();
		enter_subsection("Advection");
		{
			enter_subsection("SEDG");
			{
				declare_entry("Output directory", "/tmp",
						dealii::Patterns::DirectoryName(),
						"The name of the directory to which the output is written.");
				declare_entry("Output checkpoint interval", "1000",
						dealii::Patterns::Integer(1),
						"Write out checkpoint files every ... step.");
				declare_entry("Output solution interval", "1000",
						dealii::Patterns::Integer(1),
						"Write out solution every ... step.");
				declare_entry("Command line verbosity", "Basic",
						dealii::Patterns::Selection("Error|Basic|Full"),
						"The amount of command line output.");

			}
			leave_subsection();
			enter_subsection("5-stage Runge-Kutta");
			{

			}
			leave_subsection();
		}
		leave_subsection();
		enter_subsection("Collision");
		{

		}
		leave_subsection();

		enter_subsection("Intitialization");
		{
			enter_subsection("Iterative Initialization stop condition");
			{

			}
			leave_subsection();
		}
		leave_subsection();

		enter_subsection("Stop condition");
		{

		}
		leave_subsection();

		enter_subsection("Output");
		{

		}
		leave_subsection();
		// TODO read configuration from file
		// TODO custom configurations
		m_advectionOperatorType = Advection_SEDGMinLee;
		m_collisionType = Collision_BGKTransformed;
		m_stencilType = Stencil_D2Q9;
		m_timeIntegratorType = Integrator_RungeKutta5LowStorage;
		m_fluxType = Flux_LaxFriedrichs;
		m_timeStep = 0.1;
		m_orderOfFiniteElement = 2;
		m_numberOfTimeSteps = 100;
		m_dQScaling = 1.0;
		m_outputDirectory = "/tmp/natrium";
		setOutputFlags(out_CommandLineBasic | out_VectorFields);
		m_restart = false;
		m_distributionInitType = Equilibrium;
		m_maxDistributionInitIterations = 10000;
		m_stopDistributionInitResidual = 1e-6;
		m_outputVectorFieldsEvery = 10;
		m_outputCheckpointEvery = 500;
	}
	;

	/// destructor
	virtual ~SolverConfiguration() {
	}
	;

	/**
	 * @short Check if the configuration is consistent
	 */
	void checkConfiguration() {
		/// Test writing permission on output directory
		/// If no restart: Check if checkpoint exists. If yes -> ask for overwrite

	}

	/**
	 * @short Check if the problem definition is in accordance with the solver configuration
	 *
	 * @param[in] cFDProblem Shared pointer to a problem description
	 *
	 * @throws ... //TODO implement custom exception
	 */
	void checkProblem(shared_ptr<ProblemDescription<2> > cFDProblem) {
		//TODO: implement the checkProblem function
	}

	/**
	 * @short Check if the problem definition is in accordance with the solver configuration
	 *
	 * @param[in] cFDProblem Shared pointer to a problem description
	 *
	 * @throws ... //TODO implement custom exception
	 */
	void checkProblem(shared_ptr<ProblemDescription<3> > cFDProblem) {
		//TODO: implement the checkProblem function
	}

	CollisionType getCollisionType() const {
		return m_collisionType;
	}

	void setCollisionType(CollisionType collisionType) {
		m_collisionType = collisionType;
	}

	StencilType getStencilType() const {
		return m_stencilType;
	}

	void setStencilType(StencilType stencilType) {
		m_stencilType = stencilType;
	}

	double getTimeStep() const {
		return m_timeStep;
	}

	void setTimeStep(double timeStep) {
		m_timeStep = timeStep;
	}

	size_t getOrderOfFiniteElement() const {
		return m_orderOfFiniteElement;
	}

	void setOrderOfFiniteElement(size_t orderOfFiniteElement) {
		m_orderOfFiniteElement = orderOfFiniteElement;
	}

	TimeIntegratorType getTimeIntegratorType() const {
		return m_timeIntegratorType;
	}

	void setTimeIntegratorType(TimeIntegratorType timeIntegratorType) {
		m_timeIntegratorType = timeIntegratorType;
	}

	AdvectionOperatorType getAdvectionOperatorType() const {
		return m_advectionOperatorType;
	}

	void setAdvectionOperatorType(AdvectionOperatorType advectionOperatorType) {
		m_advectionOperatorType = advectionOperatorType;
	}

	size_t getNumberOfTimeSteps() const {
		return m_numberOfTimeSteps;
	}

	void setNumberOfTimeSteps(size_t numberOfTimeSteps) {
		m_numberOfTimeSteps = numberOfTimeSteps;
	}

	FluxType getFluxType() const {
		return m_fluxType;
	}

	void setFluxType(FluxType fluxType) {
		m_fluxType = fluxType;
	}

	const std::string& getOutputDirectory() const {
		return m_outputDirectory;
	}

	void setOutputDirectory(const std::string& outputDirectory) {
		// TODO create directory; check mod
		this->m_outputDirectory = outputDirectory;
	}

	double getDQScaling() const {
		return m_dQScaling;
	}

	void setDQScaling(double dQScaling) {
		m_dQScaling = dQScaling;
	}

	int getOutputFlags() const {
		return m_outputFlags;
	}

	void setOutputFlags(int outputFlags) {
		m_outputFlags = outputFlags;
		// if Complete log, then switch on all commandline flags
		if ((out_CommandLineFull & m_outputFlags) != 0) {
			m_outputFlags |= out_CommandLineBasic;
			m_outputFlags |= out_CommandLineError;
		}
		// if base, then switch on errors
		if ((out_CommandLineBasic & m_outputFlags) != 0) {
			m_outputFlags |= out_CommandLineError;
		}
		// redefine Logging stream
		/*std::stringstream logFile;
		 if ((out_LogFile & m_outputFlags) != 0){
		 logFile << getOutputDirectory() << "/natrium.log";
		 } else {
		 logFile << "";
		 }
		 if ((out_CommandLineFull & m_outputFlags) != 0) {
		 Logging::FULL = Logging::makeTeeStream(true, false, false, logFile.str());
		 } else {
		 Logging::FULL = Logging::makeTeeStream(false, false, false, logFile.str());
		 }
		 if ((out_CommandLineBasic & m_outputFlags) != 0) {
		 Logging::BASIC = Logging::makeTeeStream(true, false, false, logFile.str());
		 } else {
		 Logging::BASIC = Logging::makeTeeStream(false, false, false, logFile.str());
		 }
		 if ((out_CommandLineError & m_outputFlags) != 0) {
		 Logging::ERROR = Logging::makeTeeStream(true, false, false, logFile.str());
		 } else {
		 Logging::ERROR = Logging::makeTeeStream(false, false, false, logFile.str());
		 }*/

	}

	bool isRestart() const {
		return m_restart;
	}

	void setRestart(bool restart) {
		m_restart = restart;
	}

	DistributionInitType getDistributionInitType() const {
		return m_distributionInitType;
	}

	void setDistributionInitType(DistributionInitType distributionInitType) {
		m_distributionInitType = distributionInitType;
	}

	size_t getMaxDistributionInitIterations() const {
		return m_maxDistributionInitIterations;
	}

	void setMaxDistributionInitIterations(
			size_t maxDistributionInitIterations) {
		m_maxDistributionInitIterations = maxDistributionInitIterations;
	}

	double getStopDistributionInitResidual() const {
		return m_stopDistributionInitResidual;
	}

	void setStopDistributionInitResidual(double stopDistributionInitResidual) {
		m_stopDistributionInitResidual = stopDistributionInitResidual;
	}

	size_t getOutputCheckpointEvery() const {
		return m_outputCheckpointEvery;
	}

	void setOutputCheckpointEvery(size_t outputCheckpointEvery) {
		m_outputCheckpointEvery = outputCheckpointEvery;
	}

	size_t getOutputVectorFieldsEvery() const {
		return m_outputVectorFieldsEvery;
	}

	void setOutputVectorFieldsEvery(size_t outputVectorFieldsEvery) {
		m_outputVectorFieldsEvery = outputVectorFieldsEvery;
	}
};

} /* namespace natrium */
#endif /* SOLVERCONFIGURATION_H_ */
