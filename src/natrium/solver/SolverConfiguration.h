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

//////////////////////////////
// DECLARE ALL SELECTIONS ////
//////////////////////////////

/**
 * @short Implemented streaming data types
 */
enum AdvectionSchemeName {
	SEDG
};

/**
 * @short Implemented collision models
 */
enum CollisionSchemeName {
	BGK_WITH_TRANSFORMED_DISTRIBUTION_FUNCTIONS // Collision for the transformed distribution function as defined in MinLee2011
};

// StencilType defined in BoltzmannModel.h

/**
 * @short Implemented time integrators
 */
enum TimeIntegratorName {
	RUNGE_KUTTA_5STAGE
};

/**
 * @short the numerical flux used to calculate the advection operator
 */
enum FluxTypeName {
	LAX_FRIEDRICHS, CENTRAL
};

/**
 * @short the initialization procedure for the distribution functions
 */
enum InitializationSchemeName {
	EQUILIBRIUM, // Distribute with equilibrium functions
	ITERATIVE // Distribute with iterative procedure; enforces consistent initial conditions
};

//////////////////////////////
// EXCEPTION CLASS        ////
//////////////////////////////

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
	ConfigurationException(const string& msg) :
			message(msg) {
	}
	~ConfigurationException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

//////////////////////////////
// MAIN CLASS             ////
//////////////////////////////

/** @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 *  @note The class is a subclass of dealii::ParameterHandler and functions as a wrapper.
 */
class SolverConfiguration: public dealii::ParameterHandler {
private:

	//////////////////////////////
	// MEMBER VARIABLES       ////
	//////////////////////////////

	// SECTION ADVECTION
	/// Streaming data type (e.g. MinLee2011)
	AdvectionSchemeName m_advectionScheme;
	/// Time Integrator type (e.g. RK5LowStorage)
	TimeIntegratorName m_timeIntegrator;

	// --- SUBSECTION SEDG
	/// Numerical Flux
	FluxTypeName m_SEDGFluxType;
	/// Order of Finite Element
	size_t m_SEDGOrderOfFiniteElement;

	// SECTION COLLISION
	bool m_collisionOnBoundaryNodes;
	/// Collision type (e.g. BGKTransformed)
	CollisionSchemeName m_collisionScheme;

	// SECTION GENERAL
	/// Stencil type (e.g. D2Q9)
	StencilType m_stencilType;
	/// scaling of the difference stencil for the discrete particle velocity
	double m_stencilScaling;
	/// Switch output off?
	bool m_switchOutputOff;
	/// Time step size
	double m_timeStepSize;

	// SECTION INITIALIZATION
	InitializationSchemeName m_initializationScheme;
	bool m_restartAtLastCheckpoint;

	// --- SUBSECTION ITERATIVE INITIALIZATION STOP CONDITION
	/// max initialization iterations
	size_t m_iterativeInitializationNumberOfIterations;
	/// stop condition for the density residual of the initialization procedure
	double m_iterativeInitializationResidual;

	// SECTION OUTPUT
	/// degree of output on the command line
	size_t m_commandLineVerbosity;
	/// output frequency checkpoints
	size_t m_outputCheckpointInterval;
	/// Output directory
	std::string m_outputDirectory;
	/// Output frequency solution vector
	size_t m_outputSolutionInterval;
	/// Indicates whether to write a log file
	bool m_writeALogFile;

	// SECTION STOP CONDITION
	/// Number of time steps
	size_t m_numberOfTimeSteps;

public:

	/**
	 * @short Constructor -- initializes parameters with their default values
	 */
	SolverConfiguration();

	/**
	 * @short Constructor -- initializes parameters from an xmlfile
	 * @throws ConfigurationException if the XML file has wrong format
	 */
	SolverConfiguration(const std::string& XMLfilename);

	/// destructor
	virtual ~SolverConfiguration() {
	}
	;

	/**
	 * @short wrapper function for ParameterHandler::read_input; directing cerr into a C++-Exception
	 **/
	void readFromTextFile(const std::string & filename, const bool optional =
			false, const bool write_stripped_file = false);

	/**
	 * @short wrapper function for ParameterHandler::read_input_from_xml; directing cerr into a C++-Exception
	 **/
	void readFromXMLFile(const std::string & filename);

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

	/*void setOutputFlags(int outputFlags) {
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
	 std::stringstream logFile;
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
	 }

	 }*/

	//////////////////////////////////
	// GETTER AND SETTER -------  ////
	// WRAPPED AROUND THE DEAL.II ////
	// PARMETER HANDLER CLASS     ////
	//////////////////////////////////
	AdvectionSchemeName getAdvectionScheme() {
		enter_subsection("Advection");
		string advectionScheme = get("Advection scheme");
		leave_subsection();
		if ("SEDG" == advectionScheme) {
			return SEDG;
		} else {
			std::stringstream msg;
			msg << "Unknown advection scheme '" << advectionScheme
					<< " '. The implementation of AdvectionSchemeName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setAdvectionScheme(AdvectionSchemeName advectionScheme) {
		enter_subsection("Advection");
		switch (advectionScheme) {
		case SEDG: {
			set("Advection scheme", "SEDG");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown advection scheme; index. " << advectionScheme
					<< " in enum AdvectionScheme. The constructor of SolverConfiguration might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
	}

	bool isCollisionOnBoundaryNodes() {
		enter_subsection("Collision");
		bool collisionOnBoundaryNodes;
		try {
			collisionOnBoundaryNodes = get_bool("Collision on boundary nodes");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Collision on boundary nodes' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return collisionOnBoundaryNodes;
	}

	void setCollisionOnBoundaryNodes(bool collisionOnBoundaryNodes) {
		m_collisionOnBoundaryNodes = collisionOnBoundaryNodes;
	}

	CollisionSchemeName getCollisionScheme() const {
		return m_collisionScheme;
	}

	void setCollisionScheme(CollisionSchemeName collisionScheme) {
		m_collisionScheme = collisionScheme;
	}

	size_t getCommandLineVerbosity() const {
		return m_commandLineVerbosity;
	}

	void setCommandLineVerbosity(size_t commandLineVerbosity) {
		m_commandLineVerbosity = commandLineVerbosity;
	}

	InitializationSchemeName getInitializationScheme() const {
		return m_initializationScheme;
	}

	void setInitializationScheme(
			InitializationSchemeName initializationScheme) {
		m_initializationScheme = initializationScheme;
	}

	size_t getIterativeInitializationNumberOfIterations() const {
		return m_iterativeInitializationNumberOfIterations;
	}

	void setIterativeInitializationNumberOfIterations(
			size_t iterativeInitializationNumberOfIterations) {
		m_iterativeInitializationNumberOfIterations =
				iterativeInitializationNumberOfIterations;
	}

	double getIterativeInitializationResidual() const {
		return m_iterativeInitializationResidual;
	}

	void setIterativeInitializationResidual(
			double iterativeInitializationResidual) {
		m_iterativeInitializationResidual = iterativeInitializationResidual;
	}

	size_t getNumberOfTimeSteps() const {
		return m_numberOfTimeSteps;
	}

	void setNumberOfTimeSteps(size_t numberOfTimeSteps) {
		m_numberOfTimeSteps = numberOfTimeSteps;
	}

	size_t getOutputCheckpointInterval() const {
		return m_outputCheckpointInterval;
	}

	void setOutputCheckpointInterval(size_t outputCheckpointInterval) {
		m_outputCheckpointInterval = outputCheckpointInterval;
	}

	const std::string& getOutputDirectory() const {
		return m_outputDirectory;
	}

	void setOutputDirectory(const std::string& outputDirectory) {
		m_outputDirectory = outputDirectory;
	}

	size_t getOutputSolutionInterval() const {
		return m_outputSolutionInterval;
	}

	void setOutputSolutionInterval(size_t outputSolutionInterval) {
		m_outputSolutionInterval = outputSolutionInterval;
	}

	bool isRestartAtLastCheckpoint() const {
		return m_restartAtLastCheckpoint;
	}

	void setRestartAtLastCheckpoint(bool restartAtLastCheckpoint) {
		m_restartAtLastCheckpoint = restartAtLastCheckpoint;
	}

	FluxTypeName getSedgFluxType() const {
		return m_SEDGFluxType;
	}

	void setSedgFluxType(FluxTypeName sedgFluxType) {
		m_SEDGFluxType = sedgFluxType;
	}

	size_t getSedgOrderOfFiniteElement() const {
		return m_SEDGOrderOfFiniteElement;
	}

	void setSedgOrderOfFiniteElement(size_t sedgOrderOfFiniteElement) {
		m_SEDGOrderOfFiniteElement = sedgOrderOfFiniteElement;
	}

	double getStencilScaling() const {
		return m_stencilScaling;
	}

	void setStencilScaling(double stencilScaling) {
		m_stencilScaling = stencilScaling;
	}

	StencilType getStencilType() const {
		return m_stencilType;
	}

	void setStencilType(StencilType stencilType) {
		m_stencilType = stencilType;
	}

	TimeIntegratorName getTimeIntegrator() const {
		return m_timeIntegrator;
	}

	void setTimeIntegrator(TimeIntegratorName timeIntegrator) {
		m_timeIntegrator = timeIntegrator;
	}

	double getTimeStepSize() const {
		return m_timeStepSize;
	}

	void setTimeStepSize(double timeStepSize) {
		m_timeStepSize = timeStepSize;
	}

	bool isWriteALogFile() const {
		return m_writeALogFile;
	}

	void setWriteALogFile(bool writeALogFile) {
		m_writeALogFile = writeALogFile;
	}

	bool isSwitchOutputOff() const {
		return m_switchOutputOff;
	}

	void setSwitchOutputOff(bool switchOutputOff) {
		this->m_switchOutputOff = switchOutputOff;
	}
};

} /* namespace natrium */
#endif /* SOLVERCONFIGURATION_H_ */
