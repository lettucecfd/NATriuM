/**
 * @file CFDSolver.h
 * @short Central class of the CFD Simulation based on the Discrete Boltzmann Equation (DBE)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVER_H_
#define CFDSOLVER_H_

#include <exception>

#include "../problemdescription/ProblemDescription.h"

#include "../advection/AdvectionOperator.h"
#include "../advection/SEDGMinLee.h"

#include "../boltzmannmodels/BoltzmannModel.h"
#include "../boltzmannmodels/D2Q9IncompressibleModel.h"

#include "../collisionmodels/CollisionModel.h"
#include "../collisionmodels/BGKTransformed.h"

#include "../timeintegration/TimeIntegrator.h"
#include "../timeintegration/RungeKutta5LowStorage.h"

#include "../utilities/BasicNames.h"

#include "SolverConfiguration.h"

namespace natrium {

/**
 * @short Exception class for CFDSolver
 */
class CFDSolverException: public std::exception {
private:
	std::string message;
public:
	CFDSolverException(const char *msg) :
			message(msg) {
	}
	~CFDSolverException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/** @short The central class for the CFD simulation based on the DBE.
 *  @note  The CFDSolver itself is quite static but it contains interchangeable modules, e.g. for the
 *         Boltzmann model or the time integrator. By these means, a variety of different simulation
 *         methods can be covered.
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class CFDSolver {

private:

	/// particle distribution functions
	vector<distributed_vector> m_f;

	/// macroscopic density
	distributed_vector m_density;

	/// macroscopic velocity
	vector<distributed_vector> m_velocity;

	/// description of the CFD problem (boundraries, initial values, etc.)
	shared_ptr<ProblemDescription<dim> > m_problemDescription;

	/// global streaming data
	shared_ptr<AdvectionOperator<dim> > m_advectionOperator;

	/// DdQq Boltzmann model (e.g. D2Q9)
	shared_ptr<BoltzmannModel> m_boltzmannModel;

	/// Description of the collision algorithm
	shared_ptr<CollisionModel> m_collisionModel;

	/// Time Integrator for the solution of the ODE, which stems from the space discretization
	shared_ptr<TimeIntegrator> m_timeIntegrator;

	/// Configuration of the solver
	shared_ptr<SolverConfiguration> m_configuration;

public:

	/// constructor
	/// @note: has to be inlined, if the template parameter is not made explicit
	CFDSolver(shared_ptr<SolverConfiguration> configuration,
			shared_ptr<ProblemDescription<dim> > problemDescription);

/// destructor
	virtual ~CFDSolver() {
	}
	;

	/**
	 * @short Advection in all directions
	 */
	void stream();

	/**
	 *  @short Low-level collide function
	 */
	void collide();

	/**
	 * @short reassembly of all matrices
	 */
	void reassemble();

	/**
	 * @short run CFD solver
	 */
	void run();

	/**
	 * @short create output data and write to file
	 */
	void output(size_t iteration);

	const distributed_vector& getDensity() const {
		return m_density;
	}

	const vector<distributed_vector>& getF() const {
		return m_f;
	}

	const vector<distributed_vector>& getVelocity() const {
		return m_velocity;
	}
}
;

} /* namespace natrium */

#endif /* CFDSOLVER_H_ */
