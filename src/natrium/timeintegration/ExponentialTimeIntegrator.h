/**
 * @file ExponentialTimeIntegrator.h
 * @short Exponential time integration scheme for the solution of f' = L*f.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef EXPONENTIALTIMEINTEGRATOR_H_
#define EXPONENTIALTIMEINTEGRATOR_H_

#include "TimeIntegrator.h"

namespace natrium {

/** @short Exponential time integration scheme for the solution of f' = L*f, as used in
 *         Uga etal. (2012) Spectral-element discontinuous Galerkin lattice Boltzmann simulation
 *         of flow past two cylinders in tandem with an exponential time integrator, CMWA 65 pp. 239-251
 */
class ExponentialTimeIntegrator: TimeIntegrator {
public:

	/// constructor
	ExponentialTimeIntegrator();

	/// destructor
	virtual ~ExponentialTimeIntegrator();
};

} /* namespace natrium */
#endif /* EXPONENTIALTIMEINTEGRATOR_H_ */
