/**
 * @file TimeIntegrator.h
 * @short Abstract class for time integrationof ordinary differential equations (ODEs).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TIMEINTEGRATOR_H_
#define TIMEINTEGRATOR_H_

namespace natrium {

/** @short Abstract class for time integration (solution of ordinary differential equations (ODE)).
 *  @note  The ODEs arise from the space discretization of the DBE, which is basically a partial differential equation (PDE).
 *         By application of a space discretization scheme (like discontinuous Galerkin or standard FEM methods)
 *         the space derivatives in the PDE are replaced with arithmetic expressions.
 *         The only remaining derivative is then the time derivative, which makes the equation an ODE.
 *         The latter can be solved using classical time integration methods like Runge-Kutta or Adams-Moulton.
 */
class TimeIntegrator {
public:

	/// constructor
	TimeIntegrator();

	/// destructor
	virtual ~TimeIntegrator();
};

} /* namespace natrium */
#endif /* TIMEINTEGRATOR_H_ */
