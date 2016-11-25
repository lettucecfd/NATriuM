/**
 * @file PoiseuilleFlow3D.h
 * @short Channel flow in 2D
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PoiseuilleFlow3D_H_
#define PoiseuilleFlow3D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "../problemdescription/Benchmark.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Channel Flow
 *  The domain is [0,5]x[0,1].
 */
class PoiseuilleFlow3D: public Benchmark<3> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class AnalyticVelocity: public dealii::Function<3> {
	private:
		PoiseuilleFlow3D* m_flow;
	public:
		AnalyticVelocity(PoiseuilleFlow3D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;
	};

	/// constructor
	PoiseuilleFlow3D(double viscosity, size_t refinementLevel, double u_bulk = 1.0,
			double height = 1.0, double width = 2.0, double length = 2.0, bool is_periodic = true);

	/// destructor
	virtual ~PoiseuilleFlow3D();


	virtual double getCharacteristicVelocity() const {
		return m_uBulk;
	}

	virtual void refine(Mesh<3>& mesh){
		// Refine grid
		mesh.refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<3>& ){

	}
	virtual bool isCartesian(){
		return true;
	}
private:

	double m_uBulk;

	double m_uMax;

	size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<3> > makeGrid(double height, double width, double length);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<3> > makeBoundaries(bool is_periodic);



};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
