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
enum StreamingDataType {
	Streaming_MinLee2011
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



/** @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 *  @tparam dim The dimension of the flow (2 or 3).
 */
class SolverConfiguration {
private:

	/// Streaming data type (e.g. MinLee2011)
	StreamingDataType m_streamingDataType;

	/// Collision type (e.g. BGKTransformed)
	CollisionType m_collisionType;

	/// Stencil type (e.g. D2Q9)
	StencilType m_stencilType;

	/// Time Integrator type (e.g. RK5LowStorage)
	TimeIntegratorType m_timeIntegratorType;

	/// Time step size
	double m_timeStep;

	/// Order of finite element
	size_t m_orderOfFiniteElement;

	///

public:

	/// constructor
	SolverConfiguration(){
		// TODO read configuration from file
		// TODO custom configurations
		m_streamingDataType = Streaming_MinLee2011;
		m_collisionType = Collision_BGKTransformed;
		m_stencilType = Stencil_D2Q9;
		m_timeIntegratorType = Integrator_RungeKutta5LowStorage;
		m_timeStep = 1.0;
		m_orderOfFiniteElement = 1;
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

	StreamingDataType getStreamingDataType() const {
		return m_streamingDataType;
	}

	void setStreamingDataType(StreamingDataType streamingDataType) {
		m_streamingDataType = streamingDataType;
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
};

} /* namespace natrium */
#endif /* SOLVERCONFIGURATION_H_ */
