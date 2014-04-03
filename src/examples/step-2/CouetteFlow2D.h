/**
 * @file CouetteFlow2D.h
 * @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COUETTEFLOW2D_H_
#define COUETTEFLOW2D_H_

#include "deal.II/grid/tria.h"

#include "problemdescription/ProblemDescription.h"
#include "utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

/** @short Description of a simple Couette Flow (regular channel flow in square domain).
 *  The domain is [0,1]^2. The top plate is moved with constant velocity. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 *  @note The analytic solution is obtained by a formula stated in Min and Lee (2011).
 */
class CouetteFlow2D: public ProblemDescription<2> {
public:

	/// constructor
	CouetteFlow2D(double viscosity, double topPlateVelocity, size_t refinementLevel, double L=1.0);

	/// destructor
	virtual ~CouetteFlow2D();

	/*
	 * @short analytic solution of the Taylor-Green vortex, first component of velocity vector
	 */
	double analyticVelocity1(const dealii::Point<2>& x, double t) const ;

	/**
	 * @short analytic solution of the Taylor-Green vortex, second component of velocity vector
	 */
	double analyticVelocity2(const dealii::Point<2>& x, double t) const ;

	void applyInitialDensities(
			distributed_vector& initialDensities,
			vector<dealii::Point<2> >& supportPoints) const {
		for (size_t i = 0; i < initialDensities.size(); i++) {
			initialDensities(i) = 1.0;
		}
	}

	void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			vector<dealii::Point<2> >& supportPoints) const {
		assert(
				initialVelocities.at(0).size()
						== initialVelocities.at(1).size());
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
			initialVelocities.at(0)(i) = 0.0;
			initialVelocities.at(1)(i) = 0.0;
		}
	}

private:

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > makeGrid(double L, size_t refinementLevel);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries(double topPlateVelocity);

};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
