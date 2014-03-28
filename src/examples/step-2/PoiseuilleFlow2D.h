/**
 * @file PoiseuilleFlow2D.h
 * @short Channel flow in 2D
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TaylorGreenVortex2D_H_
#define TaylorGreenVortex2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "problemdescription/ProblemDescription.h"
#include "utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

class BoundaryDensity: public dealii::Function<2> {
public:
	virtual double value(const dealii::Point<2> &p,
			const unsigned int component = 0) const {
		return 1;
	}
};
class BoundaryVelocity: public dealii::Function<2> {
private:
	double m_xVelocity;
public:
	BoundaryVelocity(double xVelocity) :
			m_xVelocity(xVelocity) {
	}
	virtual void vector_value(const dealii::Point<2> &p,
			dealii::Vector<double> &values) const {
		values(0) = m_xVelocity;
		values(1) = 0.0;
	}

	double getXVelocity() const {
		return m_xVelocity;
	}
};

/** @short Description of a simple Channel Flow
 *  The domain is [0,5]x[0,1].
 */
class PoiseuilleFlow2D: public ProblemDescription<2> {
public:

	/// constructor
	PoiseuilleFlow2D(double viscosity, size_t refinementLevel);

	/// destructor
	virtual ~PoiseuilleFlow2D();

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialDensities(distributed_vector& initialDensities,
			vector<dealii::Point<2> >& supportPoints) const;

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			vector<dealii::Point<2> >& supportPoints) const;

	/**
	 * @short analytic solution of the Taylor-Green vortex, first component of velocity vector
	 */
	double analyticVelocity1(const dealii::Point<2>& x, double t) const;

	/**
	 * @short analytic solution of the Taylor-Green vortex, second component of velocity vector
	 */
	double analyticVelocity2(const dealii::Point<2>& x, double t) const;

	virtual double getCharacteristicVelocity() const {
		return 0.1/sqrt(3);
	}

private:

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > makeGrid(size_t refinementLevel);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries();

};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
