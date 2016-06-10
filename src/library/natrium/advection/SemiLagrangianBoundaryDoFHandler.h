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

enum BoundaryHitType {
	LinearRhoU
};

template<size_t dim>
struct BoundaryHit {

	BoundaryHit(const dealii::Point<dim>& coord, double t,
			const dealii::Tensor<1, dim>& n, BoundaryHitType type,
			const typename dealii::DoFHandler<dim>::cell_iterator& c, size_t out_direction) {
		coordinates = coord;
		time = t;
		faceNormal = n;
		boundaryType = type;
		c = cell;
		outgoingDirection = out_direction;
		f_out = -1e20;
	}

	// global information
	dealii::Point<dim> coordinates;
	double time;
	dealii::Tensor<1, dim> faceNormal;
	BoundaryHitType boundaryType;
	typename dealii::DoFHandler<dim>::cell_iterator cell;

	// outgoing distribution (only one!)
	size_t outgoingDirection;
	double f_out;

	// incoming distributions
	vector<size_t> incomingDirections;
	vector<dealii::TrilinosWrappers::internal::VectorReference> incomingDistributions;

};

template<size_t dim>
class SemiLagrangianBoundaryDoFHandler {
private:

	dealii::IndexSet m_dofsThatDependOnBoundaryValues;

	dealii::TrilinosWrappers::Vector m_boundaryValues;

	std::vector<BoundaryHit<dim> > m_boundaryPoints;

	dealii::TrilinosWrappers::Vector m_secondaryBoundaryValues;

	std::vector<BoundaryHit<dim> > m_secondaryBoundaryPoints;

public:
	SemiLagrangianBoundaryDoFHandler();
	virtual ~SemiLagrangianBoundaryDoFHandler();

	void calculateBoundaryValues(){}
	void applyBoundaryValues(distributed_block_vector& f){}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_ */
