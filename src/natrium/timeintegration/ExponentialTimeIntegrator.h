#ifndef EXPONENTIALTIMEINTEGRATOR_H_
#define EXPONENTIALTIMEINTEGRATOR_H_

#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

#include "TimeIntegrator.h"
#include "../utilities/BasicNames.h"


namespace natrium {

/** @short Exponential time integration scheme for the solution of f' = L*f, as used in
 *         Uga etal. (2012) Spectral-element discontinuous Galerkin lattice Boltzmann simulation
 *         of flow past two cylinders in tandem with an exponential time integrator, CMWA 65 pp. 239-251
 */
template <class MATRIX, class VECTOR> class ExponentialTimeIntegrator : public TimeIntegrator <MATRIX, VECTOR> {

private:


public:

	/// constructor
	ExponentialTimeIntegrator(double timeStepSize, size_t problemSize);

	ExponentialTimeIntegrator(double timeStepSize, size_t problemSize, size_t numberOfBlocks);

	/// destructor
	virtual ~ExponentialTimeIntegrator(){};

	virtual void step(VECTOR& vector,
				const MATRIX& systemMatrix, const VECTOR& systemVector);


};

} /* namespace natrium */
#endif /* EXPONENTIALTIMEINTEGRATOR_H_ */
