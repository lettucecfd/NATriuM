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

#include "BoundaryHit.h"
#include "SemiLagrangianVectorReferenceTypes.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../problemdescription/Boundary.h"
#include "../utilities/BasicNames.h"
#include "../stencils/Stencil.h"

namespace natrium {

/**
 * @short
 */
template<size_t dim>
class SemiLagrangianBoundaryDoFHandler {
private:

	GeneralizedDoFVector m_allValues;

	std::vector<BoundaryHit<dim> > m_primaryBoundaryHits;

	std::vector<BoundaryHit<dim> > m_secondaryBoundaryHits;

	const Stencil<dim>& m_stencil;
	// boost::shared_ptr<BoundaryCollection<dim> >& m_boundaries;

public:
	SemiLagrangianBoundaryDoFHandler(const Stencil<dim>& stencil) :
			m_stencil(stencil) {

	}
	virtual ~SemiLagrangianBoundaryDoFHandler();

	/**
	 * @short add a Boundary hit to either the primary or secondary boundary values
	 * @param[in/out] boundary_hit the BoundaryHit instance. The present function fills in the incoming directions of boundary_hit.
	 * @return position of the boundary hit in the
	 */
	size_t addBoundaryHit(BoundaryHit<dim>& boundary_hit) {
		// get directions
		boundary_hit.boundary.makeIncomingDirections(boundary_hit, m_stencil);
		// push back to boundary hit vector
		if (boundary_hit.isPrimary) {
			m_primaryBoundaryHits.push_back(boundary_hit);
		} else {
			m_secondaryBoundaryHits.push_back(boundary_hit);
		}
	}
	BoundaryHit<dim>& getBoundaryHit(size_t id, bool primary){
		if (primary){
			return m_primaryBoundaryHits.at(id);
		} else {
			return m_secondaryBoundaryHits.at(id);
		}
	}
	void calculateBoundaryValues(double time_of_next_step) {
	}
	void applyBoundaryValues(distributed_block_vector& f) {
	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_ */
