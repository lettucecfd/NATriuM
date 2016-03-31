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

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class AnalyticVelocity: public dealii::Function<2> {
	private:
		PoiseuilleFlow2D* m_flow;
	public:
		AnalyticVelocity(PoiseuilleFlow2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	/// constructor
	PoiseuilleFlow2D(double viscosity, size_t refinementLevel, double u_bulk = 1.0,
			double height = 1.0, double length = 2.0, bool is_periodic = true);

	/// destructor
	virtual ~PoiseuilleFlow2D();


	virtual double getCharacteristicVelocity() const {
		return m_uBulk;
	}

	virtual void refine(){
		// Refine grid
		getMesh()->refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<2>& mesh){

	}

private:

	double m_uBulk;

	double m_uMax;

	size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(double height, double length);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries(bool is_periodic);



};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
