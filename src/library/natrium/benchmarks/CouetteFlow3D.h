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



namespace natrium {


/** @short Description of a simple Couette Flow (regular channel flow in square domain).
 *  The domain is [0,1]^3. The top plate is moved with constant velocity. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 *  @note The analytic solution is obtained by a formula stated in Min and Lee (2011).
 */
class CouetteFlow3D: public Benchmark<3> {
public:


	/* class to describe the x-component of the analytic solution */
	class AnalyticVelocity: public dealii::Function<3> {
	private:
		CouetteFlow3D* m_benchmark;
	public:
		AnalyticVelocity(CouetteFlow3D* couette) :
				m_benchmark(couette) {
		}
		virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;

	};

	/// constructor
	CouetteFlow3D(double viscosity, double topPlateVelocity, size_t refinementLevel, double L=1.0, double startTime=0.0, bool isUnstructured=false);

	/// destructor
	virtual ~CouetteFlow3D();


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
	shared_ptr<Mesh<3> > makeGrid(double L, size_t refinementLevel, bool isUnstructured = false);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<3> > makeBoundaries(double topPlateVelocity);

	double getStartTime() const {
		return m_startTime;
	}

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
