/**
 * @file TaylorGreenVortex2D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TaylorGreenVortex2D_H_
#define TaylorGreenVortex2D_H_

#include "deal.II/grid/tria.h"

#include "natrium/problemdescription/Benchmark.h"
#include "natrium/utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class TaylorGreenVortex2D: public Benchmark<2> {
public:

	/// constructor
	TaylorGreenVortex2D(double viscosity,
			size_t refinementLevel);

	/// destructor
	virtual ~TaylorGreenVortex2D();

	/**
	 * @short analytic solution of the Taylor-Green vortex
	 */
	virtual void getAnalyticVelocity(const dealii::Point<2>& x, double t, dealii::Point<2>& velocity) const;


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
