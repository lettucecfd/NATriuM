/**
 * @file SolverConfiguration.h
 * @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SOLVERCONFIGURATION_H_
#define SOLVERCONFIGURATION_H_

#include "../problemdescription/ProblemDescription.h"
#include "../boltzmannmodels/BoltzmannModel.h"
#include "../utilities/BasicNames.h"

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
	Collision_BGKTransformed	// Collision for the transformed distribution function as defined in MinLee2011
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
	Flux_LaxFriedrichs,
	Flux_Central
};

/**
 * Output flags
 */
enum OutputFlags{
	CommandLineErrors = 1,
	CommandLineBase = 2,
	CommandLineCompleteLog = 4,
	LogFile = 8,
	VectorFields = 16
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
class SolverConfiguration {
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

public:

	/// constructor
	SolverConfiguration(){
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
		m_outputDirectory = "../results/test";
		m_outputFlags = CommandLineBase | VectorFields;
	};

	/// destructor
	virtual ~SolverConfiguration(){};

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
		if (CommandLineCompleteLog & m_outputFlags){
			m_outputFlags |= CommandLineBase;
			m_outputFlags |= CommandLineErrors;
		}
		// if base, then switch on errors
		if (CommandLineBase & m_outputFlags){
			m_outputFlags |= CommandLineErrors;
		}
	}
};

} /* namespace natrium */
#endif /* SOLVERCONFIGURATION_H_ */
