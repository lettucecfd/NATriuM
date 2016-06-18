/*
 * BoundaryHit.h
 *
 *  Created on: 15.06.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_ADVECTION_BOUNDARYHIT_H_
#define LIBRARY_NATRIUM_ADVECTION_BOUNDARYHIT_H_

#include "deal.II/base/point.h"
#include "deal.II/base/tensor.h"
#include "deal.II/dofs/dof_handler.h"

#include "../utilities/BasicNames.h"
#include "../problemdescription/Boundary.h"
#include "SemiLagrangianVectorReferenceTypes.h"

namespace natrium {

///**
// * @short An enum that describes the type of the boundary.
// */
//enum BoundaryHitType {
//	LINEAR_RHO_U,
//	LINEAR_RHO_U_CORNER
//};

struct LagrangianPathDestination {
	size_t index;
	bool isBoundaryHit;
	bool isSecondary;
	size_t direction; // in case of a boundary: the outgoing direction
	LagrangianPathDestination(size_t i, bool is_hit, bool is_secondary, size_t alpha) :
			index(i), isBoundaryHit(is_hit), isSecondary(is_secondary), direction(alpha) {
	}
	LagrangianPathDestination(size_t i, size_t alpha) :
			index(i), isBoundaryHit(false), isSecondary(false), direction(alpha) {
	}
	LagrangianPathDestination(const LagrangianPathDestination& other):
		index(other.index), isBoundaryHit(other.isBoundaryHit), isSecondary(other.isSecondary), direction(other.direction){
	}
};

/**
 * @short A struct that defines all information of a boundary hit. (When and where did the distribution hit the boundary hit, ...)
 */
template<size_t dim>
struct BoundaryHit {

	/**
	 * @short Constructor
	 * @param coord coordinates of the point where the distribution hit the boundary
	 * @param t_shift time shift: t_hit := t + dt - t_shift denotes the point in time when the distribution hit the boundary (where t denots the global time step)
	 * @param n the outside unit normal vector at the boundary hit point
	 * @param bound reference to boundary object
	 * @param c present cell
	 * @param out_dof the generalized degree of freedom that this boundary hit will stream to
	 * @note Each distribution that is required to close the streaming step needs a separate BoundaryHit instance.
	 * @note The remaining fields incomingDirections, fIn have to filled after construction.
	 */
	BoundaryHit(const dealii::Point<dim>& coord, double t_shift,
			const dealii::Tensor<1, dim>& n, const Boundary<dim>& bound,
			typename dealii::DoFHandler<dim>::cell_iterator& c,
			const GeneralizedDoF& out_dof) :
			boundary(bound), out(out_dof) {
		coordinates = coord;
		time_shift = t_shift;
		faceNormal = n;
		cell = c;
	}

	/**
	 * @short copy constructor
	 */
	BoundaryHit(const BoundaryHit& other) :
			boundary(other.boundary), out(other.out) {
		coordinates = other.coordinates;
		time_shift = other.time_shift;
		faceNormal = other.faceNormal;
		cell = other.cell;
		in = other.in;
	}

	////////////////////////
	// global information //
	////////////////////////
	/**
	 * @short coordinates of the point where the Lagrangian path hits the boundary
	 */
	dealii::Point<dim> coordinates;

	/**
	 * @short time shift; the point in time when the Lagrangian path hits the boundary is defined as t+dt-time_shift
	 */
	double time_shift;

	/**
	 * @short outward unity normal vector at the boundary point
	 */
	dealii::Tensor<1, dim> faceNormal;

	/**
	 * @short the type of the boundary hit point
	 */
	const Boundary<dim>& boundary;

	/**
	 * @short cell that contains the boundary hit point
	 */
	typename dealii::DoFHandler<dim>::cell_iterator cell;

	///////////////////////////////////////
	// outgoing distribution (only one!) //
	///////////////////////////////////////
	/**
	 * @short if false, other boundary hits depend on the present boundary hit
	 */
	GeneralizedDoF out;

	///////////////////////////
	// incoming distribution //
	///////////////////////////
	/**
	 * @short vector of references to the values of incoming directions
	 */
	vector<GeneralizedDoF> in;

};

} /* end namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_BOUNDARYHIT_H_ */

