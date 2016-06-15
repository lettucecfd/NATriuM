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
//#include "deal.II/lac/constraint_matrix.h"
#include "deal.II/fe/fe_values.h"

#include "../utilities/BasicNames.h"

namespace natrium {

class Stencil;
template<size_t dim> struct BoundaryHit;


/**
 * @short  Abstract class for the description of boundaries.
 *         Base class of PeriodicBoundary, InflowBoundary, etc.
 * @tparam dim The dimension of the boundary is the dimension of the domain -1 (
 * 	       e.g. 2-dim meshes have 1-dim boundary)
 */
template<size_t dim> class Boundary {

public:

	/// constructor
	Boundary();

	/// destructor
	virtual ~Boundary();

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const {
		return false;
	}

	/** @short is the boundary a dirichlet boundary ?
	 */
	virtual bool isLinear() const {
		return false;
	}

	/**
	 * @short Calculates outgoing distribution from incoming distributions.
	 * @param[in/out] boundary_hit The boundary hit instance that contains all information about the boundary hit.
	 * @param[in] stencil the stencil (e.g. a D2Q9 instance)
	 * @param[in] time_of_next_step the physical time at the next time step (is required here to define time-dependent boundary conditions)
	 * @note This function is used by the semi-Lagrangian advection solver. Before it is called on a
	 *       BoundaryHit instance, the BoundaryHit instance must have the right incoming directions (usually filled
	 *       by makeIncomingDirections()) and the right references in fIn (has to be filled by hand -- by the
	 *       semi-Lagrangian advection solver).
	 */
	virtual void calculate(BoundaryHit<dim>& boundary_hit,
			const Stencil& stencil, double time_of_next_step) = 0;

	/**
	 * @short Resizes boundary_hit.incomingDirections and fills it in.
	 * @param[in/out] boundary_hit The boundary hit instance that contains all information about the boundary hit.
	 * @param[in] stencil the stencil (e.g. a D2Q9 instance)
	 * @note This function is used by the semi-Lagrangian advection solver
	 */
	virtual void makeIncomingDirections(BoundaryHit<dim>& boundary_hit,
			const Stencil& stencil) = 0;

};

template<size_t dim>
inline natrium::Boundary<dim>::Boundary() {
}

template<size_t dim>
inline natrium::Boundary<dim>::~Boundary() {
}

} /* namespace natrium */

#endif /* BOUNDARYDESCRIPTION_H_ */
