/**
 * @file TimeIntegrator.h
 * @short Abstract class for time integrationof ordinary differential equations (ODEs).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TIMEINTEGRATOR_H_
#define TIMEINTEGRATOR_H_

#include "../utilities/BasicNames.h"
#include "../problemdescription/BoundaryCollection.h"

namespace natrium {

/** @short Abstract class for time integration (solution of ordinary differential equations (ODE)).
 *  @note  The ODEs arise from the space discretization of the DBE, which is basically a partial differential equation (PDE).
 *         By application of a space discretization scheme (like discontinuous Galerkin or standard FEM methods)
 *         the space derivatives in the PDE are replaced with arithmetic expressions.
 *         The only remaining derivative is then the time derivative, which makes the equation an ODE.
 *         The latter can be solved using classical time integration methods like Runge-Kutta or Adams-Moulton.
 */
template<class MATRIX, class VECTOR> class TimeIntegrator {
private:

	/// size of the time step
	double m_timeStepSize;

	/// boundary collection - is required when nonlinear boundaries are present
	boost::shared_ptr<BoundaryCollection<2> > m_boundaries2D;
	boost::shared_ptr<BoundaryCollection<3> > m_boundaries3D;

public:

	/// constructor
	TimeIntegrator(double timeStepSize);

	/// destructor
	virtual ~TimeIntegrator() {

	}

	double getTimeStepSize() const {
		return m_timeStepSize;
	}

	void setTimeStepSize(double timeStepSize) {
		assert(timeStepSize > 0.0);
		m_timeStepSize = timeStepSize;
	}

	/**
	 * @short make one time integration step on f: \f[ \frac{df}{dt} = Af+b \f].
	 * @param[in/out] f Vector of degrees of freedom
	 * @param[in] systemMatrix Matrix A
	 * @param[in] systemVector Vector b
	 * @param[in] double t global time
	 * @param[in] double dt time step size. Required to interface deal.II's embedded RK methods
	 * @return new global time
	 * @note fully virtual method. Overloaded by subclasses.
	 */
	virtual double step(VECTOR& f, const MATRIX& systemMatrix,
			VECTOR& systemVector, double t = 0, double dt = 0) = 0;

	void setBoundaryCollection(
			boost::shared_ptr<BoundaryCollection<2> > boundaries) {
		m_boundaries2D = boundaries;
	}
	void setBoundaryCollection(
			boost::shared_ptr<BoundaryCollection<3> > boundaries) {
		m_boundaries3D = boundaries;
	}
	bool hasBoundary(){
		return ( (m_boundaries3D != NULL) or (m_boundaries2D != NULL) );
	}
	void updateSystemVector(){
		if (m_boundaries2D != NULL) {
			m_boundaries2D->updateNonlinearBoundaryValues();
		} else if (m_boundaries3D != NULL) {
			m_boundaries3D->updateNonlinearBoundaryValues();
		}
	}
};

} /* namespace natrium */
#endif /* TIMEINTEGRATOR_H_ */
