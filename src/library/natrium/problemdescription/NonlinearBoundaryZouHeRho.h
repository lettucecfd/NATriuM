/*
 * NonlinearBoundaryZouHeRho.h
 *
 *  Created on: 15.12.2015
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_

#include "NonlinearBoundary.h"
#include "BoundaryTools.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

template<size_t dim>
class NonlinearBoundaryZouHeRho: public NonlinearBoundary<dim> {
private:
	boost::shared_ptr<dealii::Function<dim> > m_boundaryPressure;
	size_t m_direction;
	dealii::Tensor<1,dim> m_outwardNormal;
	double m_sign;
public:
	/**
	 * @short
	 * @param[in] direction direction of the outward normal at the wall
	 * 	-# in 2D
	 * 		- 0) left
	 * 		- 1) right
	 * 		- 2) bottom
	 * 		- 3) top
	 * 	-# in 3D
	 * 		- 0) left
	 * 		- 1) right
	 * 		- 2) bottom
	 * 		- 3) top
	 * 		- 4) front
	 * 		- 5) back
	 */
	NonlinearBoundaryZouHeRho(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundary_pressure,
			size_t direction);
	virtual ~NonlinearBoundaryZouHeRho() {

	}

	virtual void updateNonlinearBoundaryValues() const;
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
	virtual void calculate(BoundaryHit<dim>& boundary_hit, const Stencil& stencil, double time_of_next_step,
			SemiLagrangianVectorAccess& f) const {
		const numeric_vector ea = stencil.getDirection(
				boundary_hit.out.getAlpha());
		natrium_errorexit("ZouHeRho-SL not implemented.");

	}

	/**
	 * @short Resizes boundary_hit.incomingDirections and fills it in.
	 * @param[in/out] boundary_hit The boundary hit instance that contains all information about the boundary hit.
	 * @param[in] stencil the stencil (e.g. a D2Q9 instance)
	 * @note This function is used by the semi-Lagrangian advection solver
	 */
	virtual void makeIncomingDirections(BoundaryHit<dim>& boundary_hit,
			const Stencil& stencil) const {
		assert(boundary_hit.in.size() == 0);
		natrium_errorexit("ZouHeRho-SL not implemented.");
	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_ */
