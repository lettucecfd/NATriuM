/**
 * @file ShearLayer2D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef ShearLayer2D_H_
#define ShearLayer2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class ShearLayer2D: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
		ShearLayer2D* m_flow;
	public:
		InitialVelocity(ShearLayer2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};


	/// constructor (with default cs=1/sqrt(3))
	ShearLayer2D(double viscosity,
			size_t refinement_level, double u0, double kappa);

	/// destructor
	virtual ~ShearLayer2D();

	virtual void refine(){
		// Refine grid
		getMesh()->refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<2>& ){

	}

private:

	double m_u0;
	double m_kappa;
	size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid();

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();



};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
