/**
 * @file TimeIntegrator.h
 * @short Abstract class for time integrationof ordinary differential equations (ODEs).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TIMEINTEGRATOR_H_
#define TIMEINTEGRATOR_H_

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Abstract class for time integration (solution of ordinary differential equations (ODE)).
 *  @note  The ODEs arise from the space discretization of the DBE, which is basically a partial differential equation (PDE).
 *         By application of a space discretization scheme (like discontinuous Galerkin or standard FEM methods)
 *         the space derivatives in the PDE are replaced with arithmetic expressions.
 *         The only remaining derivative is then the time derivative, which makes the equation an ODE.
 *         The latter can be solved using classical time integration methods like Runge-Kutta or Adams-Moulton.
 */
template <class MATRIX, class VECTOR> class TimeIntegrator {
private:

	/// size of the time step
	double m_timeStepSize;

public:

	/// constructor
	TimeIntegrator(double timeStepSize);

	/// destructor
	virtual ~TimeIntegrator(){

	}

	double getTimeStepSize() const {
		return m_timeStepSize;
	}

	void setTimeStepSize(double timeStepSize) {
		assert(timeStepSize > 0.0);
		m_timeStepSize = timeStepSize;
	}

	/**
	 * @short make one time integration step on vector
	 *        using the system matrix
	 */
	virtual double step(VECTOR& f, const MATRIX& systemMatrix, const VECTOR& systemVector, double t = 0, double dt = 0) = 0;
};

} /* namespace natrium */
#endif /* TIMEINTEGRATOR_H_ */
