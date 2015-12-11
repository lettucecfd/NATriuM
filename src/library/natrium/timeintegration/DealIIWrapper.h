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
#include "../utilities/BasicNames.h"
#include "../solver/SolverConfiguration.h"

namespace natrium {

/**
 * @short Wrapper class to deal.II's built-in time stepping routines
 * @short This implementation supports all time stepping schemes in deal Release 8.2.1,
 *        which includes Explicit, Implicit and Embedded explicit Runge Kutta methods.
 */
template<class MATRIX, class VECTOR>
class DealIIWrapper: public TimeIntegrator<MATRIX, VECTOR> {
private:

	MATRIX const * m_systemMatrix;
	VECTOR const * m_systemVector;

	int m_solver;
	const int m_iterations = 100000;
	const double m_tol = 1e-6;


	MATRIX m_tmpMatrix;

	/// time stepping scheme in dealii
	boost::shared_ptr<dealii::TimeStepping::RungeKutta<VECTOR> > m_dealIIRKStepper;

	/// time stepping scheme, for embedded RK methods, additional methods need to be called,
	/// which requires an extra pointer
	boost::shared_ptr<dealii::TimeStepping::EmbeddedExplicitRungeKutta<VECTOR> > m_dealIIRKEmbedded;


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

	/// constructor
	DealIIWrapper(const double timeStepSize, const DealIntegratorName rkScheme, const DealSolverName linearSolver);
	DealIIWrapper(const double timeStepSize,
		const DealIntegratorName rkScheme, const DealSolverName linearSolver,	double coarsen_param, double refine_param, double min_delta, double max_delta, double refine_tol, double coarsen_tol);
	/// destructor
	virtual ~DealIIWrapper() {
	}
	;

	/**
	 * @short make one time integration step on f: \f[ \frac{df}{dt} = Af+b \f].
	 * @param[in/out] f Vector of degrees of freedom
	 * @param[in] systemMatrix Matrix A
	 * @param[in] systemVector Vector b
	 * @param[in] double t global time
	 * @param[in] double dt time step size. Required to interface deal.II's embedded RK methods
	 * @return new global time t + dt
	 */
	virtual double step(VECTOR& vector, const MATRIX& systemMatrix, const VECTOR& systemVector, double t = 0, double dt = TimeIntegrator<MATRIX,VECTOR>::getTimeStepSize());

	const boost::shared_ptr<dealii::TimeStepping::EmbeddedExplicitRungeKutta<VECTOR> >& getDealIIEmbedded() const {
		return m_dealIIRKEmbedded;
	}

	const boost::shared_ptr<dealii::TimeStepping::RungeKutta<VECTOR> >& getDealIIRKStepper() const {
		return m_dealIIRKStepper;
	}
};

} /* namespace natrium */



#endif /* DEALIIWRAPPER_H_ */
