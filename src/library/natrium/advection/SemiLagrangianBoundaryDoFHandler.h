/*
 * SemiLagrangianBoundaryDoFHandler.h
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 *
 *
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
 * @short Main class for handling semi-Lagrangian boundaries.
 *
 * The handling of boundary conditions in the semi-Lagrangian framework is quite complex.
 * The basic problem is that the distribution functions are coupled through the
 * boundary conditions. As a result, the boundary values can depend on other distribution
 * functions and (possibly time-dependent) prescribed macroscopic values. Moreover, the
 * boundary value has to be calculated at the exact time that a distribution function hits
 * the boundary. And finally, the Lagrangian paths can hit multiple boundaries during the
 * course of one time step. To deal with these issues, NATriuM has different structs and
 * classes that contain the information at the boundary that are described in the following:
 *
 * 	- BoundaryHit: This struct represents the calculation of a single outgoing distribution
 * 	  at the boundary. It stores all information required for this calculation:
 * 	  	-# the coordinates of the boundary point
 * 	  	-# the point in time when the distribution hits the boundary
 * 	  	-# the outward face normal
 * 	  	-# a reference to the boundary instance (which implements the function calculate())
 * 	  	-# a pointer to the cell that contains the boundary point
 * 	  	-# the address of the outgoing distributions (as an instance of GeneralizedDoF)
 * 	  	-# the incoming directions
 * 	  	-# the address of the incoming distributions (as instances of FunctionDepartureValue)
 * 	  We distinguish between primary and secondary boundary hits. Secondary boundary hits
 * 	  denotes boundary hits, whose outgoing distribution serves as an incoming distribution
 * 	  to an other boundary hit. This situation occurs when a Lagrangian path hits multiple
 * 	  boundary points. While the primary boundary hits write their outgoing distribution
 * 	  directly into the vector of DoFs, the secondary boundary hits write their outgoing
 * 	  distribution into a seperate vector of secondary boundary hits. This vector is stored
 * 	  in the class SecondaryBoundaryDoFVector.
 *
 * 	- SecondaryBoundaryDoFVector: This class wraps a distributed_vector in order to provide
 * 	  access to its elements via a dealii::TrilinosWrappers::internal::VectorReference. The
 * 	  handling of the IndexSet and parallelization are a bit crude and thus hidden behind this
 * 	  class.
 *
 * 	- OutgoingDistributionValue: This struct represents the address of the outgoing distribution (where to
 * 	  write the outgoing distribution). Depending on whether we have a primary or secondary
 * 	  BoundaryHit present, the outgoing direction can either be an element of the vector of
 * 	  DoFs or an element of the vector of secondary boundary values.
 *
 * 	- IncomingDistributionValue: This struct represents the address of an incoming distribution.
 * 	  The incoming distribution can either be an element of the vector of secondary boundary
 * 	  DoFs or it has to be calculated by some incoming DoFs with the respective shape function
 * 	  values.
 *
 * 	- LagrangianPathDestination: This struct represents the destination of a Lagrangian path.
 *    This destination can either be a BoundaryHit or a degree of freedom (that represents
 *    a point in the domain). Instances of LagrangianPathDestination are also used to refer
 *    to boundary hits in the SemiLagrangianBoundaryDoFHandler.
 *
 *  - LagrangianPathTracker: This struct follows a Lagrangian path to its departure point
 *    or a boundary point. If a boundary is hit, new Lagrangian paths have to be created
 *    to track the incoming distributions at a boundary back to their origin.
 *
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
	virtual ~SemiLagrangianBoundaryDoFHandler() {

	}

	/**
	 * @short add a Boundary hit to either the primary or secondary boundary values
	 * @param[in/out] boundary_hit the BoundaryHit instance. The present function fills in the incoming directions of boundary_hit.
	 * @return position of the boundary hit in the
	 */
	LagrangianPathDestination addBoundaryHit(BoundaryHit<dim>& boundary_hit) {
		// get directions
		boundary_hit.boundary.makeIncomingDirections(boundary_hit, m_stencil);
		// push back to boundary hit vector
		if (boundary_hit.out.isUnready()){
			boundary_hit.out.setUnready(false);
			boundary_hit.out.setSecondaryBoundaryHit(true);
		}
		if ( boundary_hit.out.isSecondaryBoundaryHit() ) {
			boundary_hit.out = m_secondaryBoundaryValues.appendSecondaryBoundaryDoF();
			m_secondaryBoundaryHits.push_back(boundary_hit);
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

		const size_t n_secondary = m_secondaryBoundaryHits.size();
		const size_t n_primary = m_primaryBoundaryHits.size();

		SemiLagrangianVectorAccess f(f_in, f_out, m_secondaryBoundaryValues);

		for (size_t i = 0; i < n_secondary; i++) {
			m_secondaryBoundaryHits[i].boundary.calculate(
					m_secondaryBoundaryHits[i], m_stencil, time_of_next_step, f);
		}
		for (size_t i = 0; i < n_primary; i++) {
			m_secondaryBoundaryHits[i].boundary.calculate(
					m_secondaryBoundaryHits[i], m_stencil, time_of_next_step, f);
		}
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_ */
