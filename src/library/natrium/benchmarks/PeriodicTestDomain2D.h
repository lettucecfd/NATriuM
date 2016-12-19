/**
 * @file PeriodicTestDomain2D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PERIODICTESTDOMAIN2D_H_
#define PERIODICTESTDOMAIN2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of 4 elements.
 */
class PeriodicTestDomain2D: public ProblemDescription<2> {
public:

	/// constructor
	PeriodicTestDomain2D(size_t globalRefinementLevel);

	/// destructor
	virtual ~PeriodicTestDomain2D();

	virtual void refine(Mesh<2>& mesh){
		// Refine grid
		mesh.refine_global(m_refinementLevel);
	}
	virtual void transform(Mesh<2>& ){

	}
	virtual bool isCartesian(){
		return true;
	}
private:

	size_t m_refinementLevel;

	/**
	 * @short create triangulation for periodic domain
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
#endif /* PeriodicTestDomain2D_H_ */
