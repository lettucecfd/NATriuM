/**
 * @file BoundaryDescription.h
 * @short Abstract class for Description of a Boundary object
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BOUNDARYDESCRIPTION_H_
#define BOUNDARYDESCRIPTION_H_

#include "deal.II/grid/tria.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/lac/constraint_matrix.h"
#include "deal.II/fe/fe_values.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short  Abstract class for the description of boundaries.
 *         Base class of PeriodicBoundary, InflowBoundary, etc.
 * @tparam dim The dimension of the boundary is the dimension of the domain -1 (
 * 	       e.g. 2-dim meshes have 1-dim boundary)
 */
template<size_t dim> class BoundaryDescription {

public:

	/// constructor
	BoundaryDescription();

	/// destructor
	virtual ~BoundaryDescription();

	/**
	 * @short Apply boundaries to the degrees of freedom.
	 *        This is the central function of the boundary description classes,
	 *        which is purely virtual (=0) in this abstract class.
	 *
	 * @param doFHandler The doFHandler associated with the mesh
	 * @param constraintMatrix matrix to which constraints are stored
	 */
	virtual void applyBoundaryValues(
			const shared_ptr<dealii::DoFHandler<2> > doFHandler,
			shared_ptr<dealii::ConstraintMatrix> constraintMatrix) const = 0;
};

template<size_t dim>
inline natrium::BoundaryDescription<dim>::BoundaryDescription() {
}

template<size_t dim>
inline natrium::BoundaryDescription<dim>::~BoundaryDescription() {
}

} /* namespace natrium */

#endif /* BOUNDARYDESCRIPTION_H_ */
