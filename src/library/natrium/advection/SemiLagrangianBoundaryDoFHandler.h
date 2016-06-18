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
#include "../solver/DistributionFunctions.h"

namespace natrium {

class SemiLagrangianBoundaryHitException: public NATriuMException {
private:
	std::string message;
public:
	SemiLagrangianBoundaryHitException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	SemiLagrangianBoundaryHitException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~SemiLagrangianBoundaryHitException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short
 */
template<size_t dim>
class SemiLagrangianBoundaryDoFHandler {
private:

	SecondaryBoundaryDoFVector m_secondaryBoundaryValues;

	std::vector<BoundaryHit<dim> > m_primaryBoundaryHits;

	std::vector<BoundaryHit<dim> > m_secondaryBoundaryHits;

	const Stencil& m_stencil;
	// boost::shared_ptr<BoundaryCollection<dim> >& m_boundaries;

public:
	SemiLagrangianBoundaryDoFHandler(const dealii::IndexSet& locally_owned_dofs,
			const Stencil& stencil) :
			m_secondaryBoundaryValues(locally_owned_dofs), m_stencil(stencil) {

	}
	virtual ~SemiLagrangianBoundaryDoFHandler();

	/**
	 * @short add a Boundary hit to either the primary or secondary boundary values
	 * @param[in/out] boundary_hit the BoundaryHit instance. The present function fills in the incoming directions of boundary_hit.
	 * @return position of the boundary hit in the
	 */
	LagrangianPathDestination addBoundaryHit(BoundaryHit<dim>& boundary_hit) {
		// get directions
		boundary_hit.boundary.makeIncomingDirections(boundary_hit, m_stencil);
		// push back to boundary hit vector
		if (boundary_hit.out.isSecondaryBoundaryHit()) {
			m_secondaryBoundaryHits.push_back(boundary_hit);
			// TODO append to 2ndary vector
			return LagrangianPathDestination(m_secondaryBoundaryHits.size() - 1,
					true, true, boundary_hit.out.getAlpha());
		} else {
			m_primaryBoundaryHits.push_back(boundary_hit);
			return LagrangianPathDestination(m_primaryBoundaryHits.size() - 1,
					true, false, boundary_hit.out.getAlpha());
		}
}
BoundaryHit<dim>& getBoundaryHit(LagrangianPathDestination a) {
	if (not a.isBoundaryHit) {
		throw SemiLagrangianBoundaryHitException(
				"You tried to retrieve a boundary hit without "
						"actually referencing to a boundary destination."
						"getBoundaryHit() is only defined for Lagragian"
						"path destinations at a boundary.");
	}
	if (a.isSecondary) {
		return m_secondaryBoundaryHits.at(a.index);
	} else {
		return m_primaryBoundaryHits.at(a.index);
	}
}
void calculateAndApplyBoundaryValues(DistributionFunctions& f_out,
		const DistributionFunctions& f_in, double time_of_next_step) {
	// SemiLagrangianVectorAccess
	// for all secondary boundary hits (in right order)
	// 		calculate secondary boundary values
	// 		write to generalized dof vector (secondary-boundary part)
	// for all primary boundary hits
	//		calculate primary boundary values
	//		write to dof vector (f)

}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_ */
