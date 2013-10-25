/**
 * @file RungeKutta5LowStorage.h
 * @short Fifth-order Runge-Kutta time integration scheme with low storage consumption.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef RUNGEKUTTA5LOWSTORAGE_H_
#define RUNGEKUTTA5LOWSTORAGE_H_

#include "TimeIntegrator.h"

namespace natrium {

/** @short Implementation of the fifth-order Runge-Kutta time integration scheme with low storage consumption.
 *  @note  The scheme is described in Min and Lee (2011): A spectral-element discontinuous Galerkin lattice
 *         Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 */
class RungeKutta5LowStorage: public TimeIntegrator {
public:

	/// constructor
	RungeKutta5LowStorage();

	/// destructor
	virtual ~RungeKutta5LowStorage();
};

} /* namespace natrium */
#endif /* RUNGEKUTTA5LOWSTORAGE_H_ */
