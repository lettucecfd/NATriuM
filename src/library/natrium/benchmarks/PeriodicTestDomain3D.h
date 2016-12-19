/**
 * @file PeriodicTestDomain3D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PERIODICTESTDOMAIN3D_H_
#define PERIODICTESTDOMAIN3D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of 4 elements.
 */
class PeriodicTestDomain3D: public ProblemDescription<3> {
public:

	/// constructor
	PeriodicTestDomain3D(size_t globalRefinementLevel);

	/// destructor
	virtual ~PeriodicTestDomain3D();


private:

	size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<3> > makeGrid();

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<3> > makeBoundaries();

	virtual void refine(Mesh<3>& mesh){
		// Refine grid
		mesh.refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<3>& ){

	}
	virtual bool isCartesian(){
		return true;
	}
};

} /* namespace natrium */
#endif /* PeriodicTestDomain3D_H_ */
