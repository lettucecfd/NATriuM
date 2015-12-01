/**
 * @file SolverConfiguration.h
 * @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SOLVERCONFIGURATION_H_
#define SOLVERCONFIGURATION_H_

#include <ctime>
#include <fstream>
#include <sstream>

#include "deal.II/base/parameter_handler.h"

#include "../problemdescription/ProblemDescription.h"
#include "../stencils/Stencil.h"
#include "../utilities/BasicNames.h"
#include "../utilities/Logging.h"
#include "../utilities/NATriuMException.h"

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
	BGK_STANDARD, // Standard BGK collision Collision for the distribution function as defined in MinLee2011
	BGK_STANDARD_TRANSFORMED, // BGK collisions with transformed distributions, as used in Palabos
	BGK_STEADY_STATE // Steady state preconditioning by Guo et al. (2004)
};

// StencilType defined in Stencil.h

/**
 * @short Implemented time integrators
 * @note other refers to dealii integrators which are accessed through a wrapper class
 */
enum TimeIntegratorName {
	RUNGE_KUTTA_5STAGE, THETA_METHOD, EXPONENTIAL, OTHER
};

enum DealIntegratorName {
	FORWARD_EULER,
	RK_THIRD_ORDER,
	RK_CLASSIC_FOURTH_ORDER,
	BACKWARD_EULER,
	IMPLICIT_MIDPOINT,
	CRANK_NICOLSON,
	SDIRK_TWO_STAGES,
	HEUN_EULER,
	BOGACKI_SHAMPINE,
	DOPRI,
	FEHLBERG,
	CASH_KARP,
	NONE
};

enum DealSolverName {
	BICGSTAB, CG, FGMRES, GMRES, MINRES, QMRS, RELAXATION, RICHARDSON
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
class ConfigurationException: public NATriuMException {
private:
	std::string message;
public:
	ConfigurationException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	ConfigurationException(const string& msg) :
			NATriuMException(msg), message(msg) {
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
 *  @note To integrate new parameter into the framework, you have to declare the paramters in the constructor and write the getter and setter.
 *        The setter and getter have to navigate through the sections of the Parameter handler and catch the exceptions thrown by the Parameter handler.
 *        Finally you have to update the preprocessor/NATriuM_parameters.xml file, most easily run the UnitTests and copy results/NATriuM_parameters.xml there.
 *        The latter file is automatically created according to the declared parameters.
 *        Every selection-type parameter (e.g. Advection scheme) is implemented as an enum for better handling through the command line.
 */
class SolverConfiguration: public dealii::ParameterHandler {
private:

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
			true, const bool write_stripped_file = false);

	/**
	 * @short wrapper function for ParameterHandler::read_input_from_xml; directing cerr into a C++-Exception
	 **/
	void readFromXMLFile(const std::string & filename);

	/**
	 * @short Check if the configuration is consistent
	 */
	void isConsistent();

	/**
	 * @short prepare the Output directory
	 * @note If 'User interaction' is enabled, user input is requested in case of possible overwriting
	 * @throws SolverConfigurationError, if it was not possible
	 */
	void prepareOutputDirectory();

	/**
	 * @short Check if the problem definition is in accordance with the solver configuration
	 *
	 * @param[in] cFDProblem Shared pointer to a problem description
	 *
	 * @throws ... //TODO implement custom exception
	 */
	void checkProblem(boost::shared_ptr<ProblemDescription<2> >) {
		//TODO: implement the checkProblem function
	}

	/**
	 * @short Check if the problem definition is in accordance with the solver configuration
	 *
	 * @param[in] cFDProblem Shared pointer to a problem description
	 *
	 * @throws ... //TODO implement custom exception
	 */
	void checkProblem(boost::shared_ptr<ProblemDescription<3> >) {
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
					<< " '. Check your configuration file. If everything is alright, "
					<< "the implementation of AdvectionSchemeName might not be up-to-date.";
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
					<< " in enum AdvectionSchemeName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
	}

	CollisionSchemeName getCollisionScheme() {
		enter_subsection("Collision");
		string collisionScheme = get("Collision scheme");
		leave_subsection();
		if ("BGK standard" == collisionScheme) {
			return BGK_STANDARD;
		} else if ("BGK standard transformed" == collisionScheme) {
			return BGK_STANDARD_TRANSFORMED;
		} else if ("BGK steady state" == collisionScheme) {
			return BGK_STEADY_STATE;
		} else {
			std::stringstream msg;
			msg << "Unknown collision scheme '" << collisionScheme
					<< " '. Check your configuration file. If everything is alright, "
					<< "the implementation of CollisionSchemeName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setCollisionScheme(CollisionSchemeName collisionScheme) {
		enter_subsection("Collision");
		switch (collisionScheme) {
		case BGK_STANDARD: {
			set("Collision scheme", "BGK standard");
			break;
		}
		case BGK_STANDARD_TRANSFORMED: {
			set("Collision scheme", "BGK standard transformed");
			break;
		}
		case BGK_STEADY_STATE: {
			set("Collision scheme", "BGK steady state");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown collision scheme; index. " << collisionScheme
					<< " in enum CollisionSchemeName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
	}

	double getBGKSteadyStateGamma() {
		enter_subsection("Collision");
		enter_subsection("BGK parameters");
		double gamma;
		try {
			gamma = get_double("Steady state gamma");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Steady state gamma' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return gamma;
	}

	void setBGKSteadyStateGamma(double gamma) {
		enter_subsection("Collision");
		enter_subsection("BGK parameters");
		try {
			set("Steady state gamma", gamma);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << gamma
					<< " to Steady state gamma: " << e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
	}

	size_t getCommandLineVerbosity() {
		enter_subsection("Output");
		size_t commandLineVerbosity;
		try {
			commandLineVerbosity = get_integer("Command line verbosity");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Command line verbosity' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return commandLineVerbosity;
	}

	void setCommandLineVerbosity(long int commandLineVerbosity) {
		enter_subsection("Output");
		try {
			set("Command line verbosity", commandLineVerbosity);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << commandLineVerbosity
					<< " to command line verbosity: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	bool hasAnalyticSolution() {
		enter_subsection("General");
		bool hasAnalytic;
		try {
			hasAnalytic = get_bool("Has analytic solution?");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Has analytic solution?' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return hasAnalytic;
	}

	void setHasAnalyticSolution(bool hasAnalyticSolution) {
		enter_subsection("General");
		set("Has analytic solution?", hasAnalyticSolution);
		leave_subsection();
	}

	InitializationSchemeName getInitializationScheme() {
		enter_subsection("Initialization");
		string initializationScheme = get("Initialization scheme");
		leave_subsection();
		if ("Iterative" == initializationScheme) {
			return ITERATIVE;
		} else if ("Equilibrium" == initializationScheme) {
			return EQUILIBRIUM;
		} else {
			std::stringstream msg;
			msg << "Unknown initialization scheme '" << initializationScheme
					<< " '. Check your configuration file. If everything is alright, "
					<< "the implementation of InitilizationSchemeName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setInitializationScheme(
			InitializationSchemeName initializationScheme) {
		enter_subsection("Initialization");
		switch (initializationScheme) {
		case ITERATIVE: {
			set("Initialization scheme", "Iterative");
			break;
		}
		case EQUILIBRIUM: {
			set("Initialization scheme", "Equilibrium");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown initialization scheme scheme; index. "
					<< initializationScheme
					<< " in enum InitializationSchemeName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
	}

	size_t getIterativeInitializationNumberOfIterations() {
		enter_subsection("Initialization");
		enter_subsection("Iterative initialization stop condition");
		size_t nofIterations;
		try {
			nofIterations = get_integer("Number of iterations");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Number of iterations' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return nofIterations;
	}

	void setIterativeInitializationNumberOfIterations(
			long int iterativeInitializationNumberOfIterations) {
		enter_subsection("Initialization");
		enter_subsection("Iterative initialization stop condition");
		try {
			set("Number of iterations",
					iterativeInitializationNumberOfIterations);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value "
					<< iterativeInitializationNumberOfIterations
					<< " to Number of iterations: " << e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
	}

	double getIterativeInitializationResidual() {
		enter_subsection("Initialization");
		enter_subsection("Iterative initialization stop condition");
		double resid;
		try {
			resid = get_double("Residual");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not read parameter 'Residual' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return resid;
	}

	void setIterativeInitializationResidual(
			double iterativeInitializationResidual) {
		enter_subsection("Initialization");
		enter_subsection("Iterative initialization stop condition");
		try {
			set("Residual", iterativeInitializationResidual);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << iterativeInitializationResidual
					<< " to Iterative initialization residual: " << e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
	}

	size_t getNumberOfTimeSteps() {
		enter_subsection("Stop condition");
		size_t nofSteps;
		try {
			nofSteps = get_integer("Number of time steps");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Number of time steps' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return nofSteps;
	}

	void setNumberOfTimeSteps(long int numberOfTimeSteps) {
		enter_subsection("Stop condition");
		try {
			set("Number of time steps", numberOfTimeSteps);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << numberOfTimeSteps
					<< " to Number of time steps: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	double getSimulationEndTime() {
		enter_subsection("Stop condition");
		double end_time;
		try {
			end_time = get_double("Simulation end time");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Simulation end time' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return end_time;
	}

	void setSimulationEndTime(double end_time) {
		enter_subsection("Stop condition");
		try {
			set("Simulation end time", end_time);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << end_time
					<< " to Simulation end time: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	double getConvergenceThreshold() {
		enter_subsection("Stop condition");
		double threshold;
		try {
			threshold = get_double("Convergence threshold");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Convergence threshold' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return threshold;
	}

	void setConvergenceThreshold(double threshold) {
		enter_subsection("Stop condition");
		try {
			set("Convergence threshold", threshold);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << threshold
					<< " to Convergence threshold: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	size_t getOutputCheckpointInterval() {
		enter_subsection("Output");
		size_t checkpointInterval;
		try {
			checkpointInterval = get_integer("Output checkpoint interval");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Output checkpoint interval' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return checkpointInterval;
	}

	void setOutputCheckpointInterval(long int outputCheckpointInterval) {
		enter_subsection("Output");
		try {
			set("Output checkpoint interval", outputCheckpointInterval);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << outputCheckpointInterval
					<< " to Output checkpoint interval: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	size_t getOutputTableInterval() {
		enter_subsection("Output");
		size_t tableInterval;
		try {
			tableInterval = get_integer("Output table interval");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Output table interval' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return tableInterval;
	}

	void setOutputTableInterval(long int outputTableInterval) {
		enter_subsection("Output");
		try {
			set("Output table interval", outputTableInterval);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << outputTableInterval
					<< " to Output table interval: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	const std::string getOutputDirectory() {
		enter_subsection("Output");
		string outputDir;
		try {
			outputDir = get("Output directory");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Output directory' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return outputDir;
	}

	void setOutputDirectory(const std::string& outputDirectory) {
		enter_subsection("Output");
		try {
			set("Output directory", outputDirectory);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << outputDirectory
					<< " to Output directory: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	size_t getOutputSolutionInterval() {
		enter_subsection("Output");
		size_t solutionInterval;
		try {
			solutionInterval = get_integer("Output solution interval");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Output solution interval' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return solutionInterval;
	}

	void setOutputSolutionInterval(long int outputSolutionInterval) {
		enter_subsection("Output");
		try {
			set("Output solution interval", outputSolutionInterval);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << outputSolutionInterval
					<< " to Output solution interval: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	bool isRestartAtLastCheckpoint() {
		enter_subsection("Initialization");
		bool restartAtLastCheckpoint;
		try {
			restartAtLastCheckpoint = get_bool("Restart at last checkpoint?");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Restart at last checkpoint?' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return restartAtLastCheckpoint;
	}

	void setRestartAtLastCheckpoint(bool restartAtLastCheckpoint) {
		enter_subsection("Initialization");
		set("Restart at last checkpoint?", restartAtLastCheckpoint);
		leave_subsection();
	}

	FluxTypeName getSedgFluxType() {
		enter_subsection("Advection");
		enter_subsection("SEDG");
		string fluxType = get("Flux type");
		leave_subsection();
		leave_subsection();
		if ("Lax-Friedrichs" == fluxType) {
			return LAX_FRIEDRICHS;
		} else if ("Central" == fluxType) {
			return CENTRAL;
		} else {
			std::stringstream msg;
			msg << "Unknown Flux type '" << fluxType
					<< " '. Check your configuration file. If everything is alright, "
					<< "the implementation of FluxTypeName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setSedgFluxType(FluxTypeName sedgFluxType) {
		enter_subsection("Advection");
		enter_subsection("SEDG");
		switch (sedgFluxType) {
		case LAX_FRIEDRICHS: {
			set("Flux type", "Lax-Friedrichs");
			break;
		}
		case CENTRAL: {
			set("Flux type", "Central");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown Flux type; index. " << sedgFluxType
					<< " in enum FluxTypeName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
		leave_subsection();
	}

	size_t getSedgOrderOfFiniteElement() {
		enter_subsection("Advection");
		enter_subsection("SEDG");
		size_t orderOfFE;
		try {
			orderOfFE = get_integer("Order of finite element");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Order of finite element' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return orderOfFE;
	}

	void setSedgOrderOfFiniteElement(long int sedgOrderOfFiniteElement) {
		enter_subsection("Advection");
		enter_subsection("SEDG");
		try {
			set("Order of finite element", sedgOrderOfFiniteElement);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << sedgOrderOfFiniteElement
					<< " to Order of finite element: " << e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
	}

	double getStencilScaling() {
		enter_subsection("General");
		double stencilScaling;
		try {
			stencilScaling = get_double("Stencil scaling");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Stencil scaling' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return stencilScaling;
	}

	void setStencilScaling(double stencilScaling) {
		enter_subsection("General");
		try {
			set("Stencil scaling", stencilScaling);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << stencilScaling
					<< " to Stencil scaling: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	StencilType getStencil() {
		enter_subsection("General");
		string stencil = get("Stencil");
		leave_subsection();
		if ("D2Q9" == stencil) {
			return Stencil_D2Q9;
		} else if ("D3Q19" == stencil) {
			return Stencil_D3Q19;
		} else if ("D3Q15" == stencil) {
			return Stencil_D3Q15;
		} else if ("D3Q27" == stencil) {
			return Stencil_D3Q27;
		} else {
			std::stringstream msg;
			msg << "Unknown Stencil with index " << stencil
					<< "in enum StencilType. Check your configuration file. If everything is alright, "
					<< "the implementation of FluxTypeName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setStencil(StencilType stencil) {
		enter_subsection("General");
		switch (stencil) {
		case Stencil_D2Q9: {
			set("Stencil", "D2Q9");
			break;
		}
		case Stencil_D3Q19: {
			set("Stencil", "D3Q19");
			break;
		}
		case Stencil_D3Q15: {
			set("Stencil", "D3Q15");
			break;
		}
		case Stencil_D3Q27: {
			set("Stencil", "D3Q27");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown Stencil; index. " << stencil
					<< " in enum StencilType. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
	}

	TimeIntegratorName getTimeIntegrator() {
		enter_subsection("Advection");
		string integrator = get("Time integrator");
		leave_subsection();
		if ("Runge-Kutta 5-stage" == integrator) {
			return RUNGE_KUTTA_5STAGE;
		} else if ("Theta method" == integrator) {
			return THETA_METHOD;
		} else if ("Exponential" == integrator) {
			return EXPONENTIAL;
		} else if ("Other" == integrator) {
			return OTHER;
		} else {
			std::stringstream msg;
			msg << "Unknown Time integrator with index " << integrator
					<< "in enum TimeIntegratorName. Check your configuration file. If everything is alright, "
					<< "the implementation of TimeIntegratorName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setTimeIntegrator(TimeIntegratorName timeIntegrator) {
		enter_subsection("Advection");
		switch (timeIntegrator) {
		case RUNGE_KUTTA_5STAGE: {
			set("Time integrator", "Runge-Kutta 5-stage");
			break;
		}
		case THETA_METHOD: {
			set("Time integrator", "Theta method");
			break;
		}
		case EXPONENTIAL: {
			set("Time integrator", "Exponential");
			break;
		}
		case OTHER: {
			set("Time integrator", "Other");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown Time integrator; index. " << timeIntegrator
					<< " in enum TimeIntegratorName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
	}

	double getThetaMethodTheta() {
		enter_subsection("Advection");
		enter_subsection("Theta method");
		double theta;
		try {
			theta = get_double("Theta");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not read parameter 'Theta' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return theta;
	}

	void setThetaMethodTheta(double theta) {
		enter_subsection("Advection");
		enter_subsection("Theta method");
		try {
			set("Theta", theta);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << theta
					<< " to Theta in Theta Method: " << e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
	}

	DealIntegratorName getDealIntegrator() {
		enter_subsection("Advection");
		enter_subsection("Deal.II integrator");
		string integrator = get("Runge Kutta scheme");
		leave_subsection();
		leave_subsection();
		if ("None" == integrator) {
			return NONE;
		} else if ("Forward Euler" == integrator) {
			return FORWARD_EULER;
		} else if ("RK 3rd order" == integrator) {
			return RK_THIRD_ORDER;
		} else if ("RK Classic 4th order" == integrator) {
			return RK_CLASSIC_FOURTH_ORDER;
		} else if ("Backward Euler" == integrator) {
			return BACKWARD_EULER;
		} else if ("Implicit midpoint" == integrator) {
			return IMPLICIT_MIDPOINT;
		} else if ("Crank-Nicoloson" == integrator) {
			return CRANK_NICOLSON;
		} else if ("SDIRK 2 stages" == integrator) {
			return SDIRK_TWO_STAGES;
		} else if ("Heun-Euler" == integrator) {
			return HEUN_EULER;
		} else if ("Bogacki-Shampine" == integrator) {
			return BOGACKI_SHAMPINE;
		} else if ("Dopri" == integrator) {
			return DOPRI;
		} else if ("Fehlberg" == integrator) {
			return FEHLBERG;
		} else if ("Cash-Karp" == integrator) {
			return CASH_KARP;
		} else {
			std::stringstream msg;
			msg << "Unknown Dealii time integrator with index " << integrator
					<< "in enum TimeIntegratorName. Check your configuration file. If everything is alright, "
					<< "the implementation of DealIntegratorName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setDealIntegrator(DealIntegratorName integrator) {
		enter_subsection("Advection");
		enter_subsection("Deal.II integrator");
		switch (integrator) {
		case NONE: {
			set("Runge Kutta scheme", "None");
			break;
		}
		case FORWARD_EULER: {
			set("Runge Kutta scheme", "Forward Euler");
			break;
		}
		case RK_THIRD_ORDER: {
			set("Runge Kutta scheme", "RK 3rd order");
			break;
		}
		case RK_CLASSIC_FOURTH_ORDER: {
			set("Runge Kutta scheme", "RK Classic 4th order");
			break;
		}
		case BACKWARD_EULER: {
			set("Runge Kutta scheme", "Backward Euler");
			break;
		}
		case IMPLICIT_MIDPOINT: {
			set("Runge Kutta scheme", "Implicit midpoint");
			break;
		}
		case CRANK_NICOLSON: {
			set("Runge Kutta scheme", "Crank-Nicoloson");
			break;
		}
		case SDIRK_TWO_STAGES: {
			set("Runge Kutta scheme", "SDIRK 2 stages");
			break;
		}
		case HEUN_EULER: {
			set("Runge Kutta scheme", "Heun-Euler");
			break;
		}
		case BOGACKI_SHAMPINE: {
			set("Runge Kutta scheme", "Bogacki-Shampine");
			break;
		}
		case DOPRI: {
			set("Runge Kutta scheme", "Dopri");
			break;
		}
		case FEHLBERG: {
			set("Runge Kutta scheme", "Fehlberg");
			break;
		}
		case CASH_KARP: {
			set("Runge Kutta scheme", "Cash-Karp");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown Deal.II integrator (Runge Kutta Scheme); index. "
					<< integrator
					<< " in enum DealIntegratorName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
		leave_subsection();
	}

	DealSolverName getDealLinearSolver() {
		enter_subsection("Advection");
		enter_subsection("Deal.II linear solver");
		string solver = get("Linear solver");
		leave_subsection();
		leave_subsection();

		if ("Bicgstab" == solver) {

			return BICGSTAB;
		} else if ("Cg" == solver) {
			return CG;
		} else if ("Fgmres" == solver) {
			return FGMRES;
		} else if ("Gmres" == solver) {
			return GMRES;
		} else if ("Minres" == solver) {
			return MINRES;
		} else if ("Qmrs" == solver) {
			return QMRS;
		} else if ("Relaxation" == solver) {
			return RELAXATION;
		} else if ("Richardson" == solver) {
			return RICHARDSON;
		} else {
			std::stringstream msg;
			msg << "Unknown Dealii linear solver with index " << solver
					<< "in enum LinearSolverName. Check your configuration file. If everything is alright, "
					<< "the implementation of DealSolverName might not be up-to-date.";
			throw ConfigurationException(msg.str());
		}
	}

	void setDealLinearSolver(DealSolverName solver) {
		enter_subsection("Advection");
		enter_subsection("Deal.II linear solver");
		switch (solver) {
		case BICGSTAB: {
			set("Linear solver", "Bicgstab");
			break;
		}
		case CG: {
			set("Linear solver", "Cg");
			break;
		}
		case FGMRES: {
			set("Linear solver", "Fgmres");
			break;
		}
		case GMRES: {
			set("Linear solver", "Gmres");
			break;
		}
		case MINRES: {
			set("Linear solver", "Minres");
			break;
		}
		case QMRS: {
			set("Linear solver", "Qmrs");
			break;
		}
		case RELAXATION: {
			set("Linear solver", "Relaxation");
			break;
		}
		case RICHARDSON: {
			set("Linear solver", "Richardson");
			break;
		}
		default: {
			std::stringstream msg;
			msg << "Unknown Deal.II linear solver (Linear solver); index. "
					<< solver
					<< " in enum DealLinearSolverName. The constructor of SolverConfiguration might not be up-to-date.";
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		}
		leave_subsection();
		leave_subsection();
	}

	void setEmbeddedDealIntegratorParameters(double coarsen_param,
			double refine_param, double min_delta, double max_delta,
			double refine_tol, double coarsen_tol) {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");

		set("Coarsen parameter", coarsen_param);
		set("Refinement parameter", refine_param);
		set("Minimum time step", min_delta);
		set("Maximum time step", max_delta);
		set("Refinement tolerance", refine_tol);
		set("Coarsen tolerance", coarsen_tol);

		leave_subsection();
		leave_subsection();
	}

	double getEmbeddedDealIntegratorCoarsenParameter() {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");
		double coarsen_param;
		try {
			coarsen_param = get_double("Coarsen parameter");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Coarsen parameter' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return coarsen_param;
	}

	double getEmbeddedDealIntegratorRefinementParameter() {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");
		double refine_param;
		try {
			refine_param = get_double("Refinement parameter");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Refinement parameter' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return refine_param;
	}

	double getEmbeddedDealIntegratorMinimumTimeStep() {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");
		double min_delta;
		try {
			min_delta = get_double("Minimum time step");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Minimum time step' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return min_delta;
	}

	double getEmbeddedDealIntegratorMaximumTimeStep() {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");
		double max_delta;
		try {
			max_delta = get_double("Maximum time step");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Maximum time step' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return max_delta;
	}

	double getEmbeddedDealIntegratorRefinementTolerance() {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");
		double refine_tol = 1e-8;
		try {
			refine_tol = get_double("Refinement tolerance");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Refinement tolerance' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return refine_tol;
	}

	double getEmbeddedDealIntegratorCoarsenTolerance() {
		enter_subsection("Advection");
		enter_subsection("Embedded Parameters");
		double coarsen_tol;
		try {
			coarsen_tol = get_double("Coarsen tolerance");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Coarsen tolerance' from parameters: "
					<< e.what();
			leave_subsection();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		leave_subsection();
		return coarsen_tol;
	}

	double getTimeStepSize() {
		enter_subsection("General");
		double stepSize;
		try {
			stepSize = get_double("Time step size");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not read parameter 'Time step size' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return stepSize;
	}

	void setTimeStepSize(double timeStepSize) {
		enter_subsection("General");
		try {
			set("Time step size", timeStepSize);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Could not assign value " << timeStepSize
					<< " to Time step size: " << e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
	}

	bool isWriteALogFile() {
		enter_subsection("Output");
		bool writeLogFile;
		try {
			writeLogFile = get_bool("Write a log file?");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Write a log file?' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return writeLogFile;
	}

	void setWriteALogFile(bool writeALogFile) {
		enter_subsection("Output");
		set("Write a log file?", writeALogFile);
		leave_subsection();
	}

	bool isSwitchOutputOff() {
		enter_subsection("Output");
		bool outputOff;
		try {
			outputOff = get_bool("Switch output off?");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'Switch output off?' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return outputOff;
	}

	void setSwitchOutputOff(bool switchOutputOff) {
		enter_subsection("Output");
		set("Switch output off?", switchOutputOff);
		leave_subsection();
	}

	bool isUserInteraction() {
		enter_subsection("Output");
		bool userInteract;
		try {
			userInteract = get_bool("User interaction?");
		} catch (std::exception& e) {
			std::stringstream msg;
			msg
					<< "Could not read parameter 'User interaction?' from parameters: "
					<< e.what();
			leave_subsection();
			throw ConfigurationException(msg.str());
		}
		leave_subsection();
		return userInteract;
	}

	void setUserInteraction(bool userInteract) {
		enter_subsection("Output");
		set("User interaction?", userInteract);
		leave_subsection();
	}

}
;

} /* namespace natrium */
#endif /* SOLVERCONFIGURATION_H_ */

