/*
 * DealIIWrapper.h
 *
 *  Created on: Feb 5, 2015
 *      Author: kraemer
 */

#ifndef DEALIIWRAPPER_H_
#define DEALIIWRAPPER_H_

#include "deal.II/base/time_stepping.h"
#include "deal.II/base/std_cxx11/function.h"

#include "TimeIntegrator.h"
#include "utilities/BasicNames.h"
#include "solver/SolverConfiguration.h"

namespace natrium {

/**
 * @short not yet implemented
 */
template<class MATRIX, class VECTOR>
class DealIIWrapper: public TimeIntegrator<MATRIX, VECTOR> {
private:

	MATRIX const * m_systemMatrix;
	VECTOR const * m_systemVector;

	MATRIX m_tmpMatrix;

	/// time stepping scheme in dealii
	shared_ptr<dealii::TimeStepping::RungeKutta<VECTOR> > m_dealIIRKStepper;

	/**
	 * @short function that is needed to call evolve_one_time_step
	 */
	VECTOR evaluateF(const double t, const VECTOR & f) const;

	/**
	 * @short function that is needed to call evolve_one_time_step
	 */
	VECTOR evaluateJInverse(const double t, const double tau, const VECTOR & f) const;

	/**
	 * @short function that is needed to call evolve_one_time_step
	 */
	VECTOR evaluateIdMinusTauJInverse(const double t, const double tau,
			const VECTOR & f);

public:
	DealIIWrapper(const double timeStepSize, const DealIntegratorName rkScheme);
	virtual ~DealIIWrapper() {
	}
	;
	virtual void step(VECTOR& vector, const MATRIX& systemMatrix,
			const VECTOR& systemVector);
};

} /* namespace natrium */



#endif /* DEALIIWRAPPER_H_ */
