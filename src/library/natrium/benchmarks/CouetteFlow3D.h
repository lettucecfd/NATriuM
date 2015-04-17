/**
 * @file CouetteFlow3D.h
 * @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COUETTEFLOW3D_H_
#define COUETTEFLOW3D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/Benchmark.h"
#include "../utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

/** @short Description of a simple Couette Flow (regular channel flow in square domain).
 *  The domain is [0,1]^3. The top plate is moved with constant velocity. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 *  @note The analytic solution is obtained by a formula stated in Min and Lee (2011).
 */
class CouetteFlow3D: public Benchmark<3> {
public:

	/// constructor
	CouetteFlow3D(double viscosity, double topPlateVelocity, size_t refinementLevel, double L=1.0, double startTime=0.0, bool isUnstructured=false);

	/// destructor
	virtual ~CouetteFlow3D();

	/**
	 * @short analytic solution of the Taylor-Green vortex
	 */
	virtual void getAnalyticVelocity(const dealii::Point<3>& x, double t, dealii::Point<3>& velocity) const;


	virtual double getCharacteristicVelocity() const {
		return m_topPlateVelocity;
	}

private:

	const double m_topPlateVelocity;

	const double m_startTime;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<3> > makeGrid(double L, size_t refinementLevel, bool isUnstructured = false);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<3> > makeBoundaries(double topPlateVelocity);

	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc
	{
	  double trans(const double y) const
	  {
	    return std::tanh(2*y)/tanh(2);
	  }
	  dealii::Point<3> operator() (const dealii::Point<3> &in) const
	  {
	    return dealii::Point<3> (in(0),
	                     trans(in(1)));
	  }
	};

};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
