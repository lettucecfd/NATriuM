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
	class AnalyticVelocity: public dealii::Function<2> {
	private:
		TaylorGreenVortex2D* m_flow;
	public:
		AnalyticVelocity(TaylorGreenVortex2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
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
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	/// constructor (with default cs=1/sqrt(3))
	TaylorGreenVortex2D(double viscosity,
			size_t refinementLevel, double cs = 0.57735026919, bool init_rho_analytically = false, double L = 2*M_PI);

	/// destructor
	virtual ~TaylorGreenVortex2D();

	virtual void refine(Mesh<2>& mesh){
		// Refine grid
		mesh.refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<2>& ){

	}

	virtual bool isCartesian(){
		return true;
	}

	void setHorizontalVelocity(double u){
		m_horizontalVelocity = u;
	}

	double getHorizontalVelocity() const {
		return m_horizontalVelocity;
	}

private:
	/// speed of sound
	double m_cs;

	/// initialization with analytic pressure
	bool m_analyticInit;

	size_t m_refinementLevel;

	double m_horizontalVelocity;

	double m_L;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(double L);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();



};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
