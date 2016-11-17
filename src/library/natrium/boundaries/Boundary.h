/**
 * @file Boundary.h
 * @short Abstract class for Description of a Boundary object
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BOUNDARYDESCRIPTION_H_
#define BOUNDARYDESCRIPTION_H_

#include "deal.II/grid/tria.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_values.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short  Abstract class for the description of boundaries.
 *         Base class of PeriodicBoundary, InflowBoundary, etc.
 * @tparam dim The dimension of the boundary is the dimension of the domain -1 (
 * 	       e.g. 2-dim meshes have 1-dim boundary)
 */
template<size_t dim> class Boundary {

public:

	/// constructor
	Boundary(){};

	/// destructor
	virtual ~Boundary(){};

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const = 0;

	/** @short is the boundary a linear flux boundary as in SEDG-LBM
	 */
	virtual bool isLinearFluxBoundary() const = 0;

	virtual bool isSLBoundary() const = 0;

};



} /* namespace natrium */

#endif /* BOUNDARYDESCRIPTION_H_ */
