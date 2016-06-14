/*
 * SemiLagrangianBoundaryDoFHandler.h
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_
#define LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_

#include "deal.II/base/index_set.h"
#include "deal.II/base/tensor.h"
#include "deal.II/dofs/dof_handler.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short An enum that describes the type of the boundary.
 */
enum BoundaryHitType {
	LINEAR_RHO_U
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
	 * @param type the type of the boundary, e.g. LINEAR_RHO_U
	 * @param c present cell
	 * @param out_direction the unknown distribution that is to be calculated at the boundary
	 * @note Each distribution that is required to close the streaming step needs a separate BoundaryHit instance.
	 * @note The remaining fields incomingDirections, fIn have to filled after construction.
	 */
	BoundaryHit(const dealii::Point<dim>& coord, double t_shift,
			const dealii::Tensor<1, dim>& n, BoundaryHitType type,
			typename dealii::DoFHandler<dim>::cell_iterator& c, size_t out_direction) {
		coordinates = coord;
		time_shift = t_shift;
		faceNormal = n;
		boundaryType = type;
		cell = c;
		outgoingDirection = out_direction;
		fOut = -1e20;
	}

	/**
	 * @short copy constructor
	 */
	BoundaryHit (const BoundaryHit& other){
		coordinates = other.coordinates;
		time_shift = other.time_shift;
		faceNormal = other.faceNormal;
		boundaryType = other.boundaryType;
		cell = other.cell;
		outgoingDirection = other.outgoingDirection;
		fOut = other.fOut;
		incomingDirections = other.incomingDirections;
		fIn = other.fIn;
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
	BoundaryHitType boundaryType;

	/**
	 * @short cell that contains the boundary hit point
	 */
	typename dealii::DoFHandler<dim>::cell_iterator cell;


	///////////////////////////////////////
	// outgoing distribution (only one!) //
	///////////////////////////////////////
	/**
	 * @short index of outgoing direction
	 */
	size_t outgoingDirection;
	/**
	 * @short value of the outgoing distribution (has to be calculated in each time step)
	 */
	double fOut;

	///////////////////////////
	// incoming distribution //
	///////////////////////////
	/**
	 * @short vector that contains the indices of incoming directions
	 */
	vector<size_t> incomingDirections;

	/**
	 * @short vector of references to the values of incoming directions
	 */
	vector<dealii::TrilinosWrappers::internal::VectorReference> fIn;

};

/**
 * @short
 */
template<size_t dim>
class SemiLagrangianBoundaryDoFHandler {
private:

	dealii::IndexSet m_dofsThatDependOnBoundaryValues;

	distributed_vector m_boundaryValues;

	std::vector<BoundaryHit<dim> > m_boundaryPoints;

	distributed_vector m_secondaryBoundaryValues;

	std::vector<BoundaryHit<dim> > m_secondaryBoundaryPoints;

public:
	SemiLagrangianBoundaryDoFHandler();
	virtual ~SemiLagrangianBoundaryDoFHandler();

	void addBoundaryHit(BoundaryHit<dim>& boundary_hit){

	}
	void calculateBoundaryValues(double time_of_next_step){}
	void applyBoundaryValues(distributed_block_vector& f){}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_ */
