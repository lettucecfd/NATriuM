/**
 * @file TaylorGreenVortex2D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TaylorGreenVortex2D_H_
#define TaylorGreenVortex2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/Benchmark.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class TaylorGreenVortex2D: public Benchmark<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 */
	class AnalyticVelocityU: public dealii::Function<2> {
	private:
		TaylorGreenVortex2D* m_flow;
	public:
		AnalyticVelocityU(TaylorGreenVortex2D* flow) :
				m_flow(flow) {
		}
		double value(const dealii::Point<2>& x) const;
	};

	/**
	 * @short class to describe the y-component of the analytic solution
	 */
	class AnalyticVelocityV: public dealii::Function<2> {
	private:
		TaylorGreenVortex2D* m_flow;
	public:
		AnalyticVelocityV(TaylorGreenVortex2D* flow) :
				m_flow(flow) {
		}
		double value(const dealii::Point<2>& x) const;
	};

	/**
	 * @short class to describe the y-component of the analytic solution
	 */
	class AnalyticDensity: public dealii::Function<2> {
	private:
		TaylorGreenVortex2D* m_flow;
	public:
		AnalyticDensity(TaylorGreenVortex2D* flow) :
				m_flow(flow) {
		}
		double value(const dealii::Point<2>& x) const;
	};

	/// constructor (with default cs=1/sqrt(3))
	TaylorGreenVortex2D(double viscosity,
			size_t refinementLevel, double cs = 0.57735026919, bool init_rho_analytically = false);

	/// destructor
	virtual ~TaylorGreenVortex2D();


private:
	/// speed of sound
	double m_cs;

	/// initialization with analytic pressure
	bool m_analyticInit;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Mesh<2> > makeGrid(size_t refinementLevel);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries();

};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
