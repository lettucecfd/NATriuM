/**
 * @file PoiseuilleFlow2D.h
 * @short Channel flow in 2D
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef POISEUILLEFLOW2D_H_
#define POISEUILLEFLOW2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "../problemdescription/Benchmark.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Channel Flow
 *  The domain is [0,5]x[0,1].
 */
class PoiseuilleFlow2D: public Benchmark<2> {
public:

	/// constructor
	PoiseuilleFlow2D(double viscosity, size_t refinementLevel, double u_bulk = 1.0,
			double height = 1.0, double length = 2.0, bool is_periodic = true);

	/// destructor
	virtual ~PoiseuilleFlow2D();

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<2> >& supportPoints) const;

	/**
	 * @short analytic solution of the Taylor-Green vortex
	 */
	virtual void getAnalyticVelocity(const dealii::Point<2>& x, double t,
			dealii::Point<2>& velocity) const;

	/**
	 * @short get Analytic density at one point in space and time
	 * @note Lattice Boltzmann uses an ideal EOS to couple density and pressure. In that sense, the analytic
	 * density here is rather a pressure, and as such calculated from the analytic solution for p.
	 */
	virtual double getAnalyticDensity(const dealii::Point<2>& x,
			double t) const;

	virtual double getCharacteristicVelocity() const {
		return m_uBulk;
	}

private:

	double m_uBulk;

	double m_uMax;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Mesh<2> > makeGrid(size_t refinementLevel,
			double height, double length, bool is_periodic);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries(bool is_periodic);

};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
