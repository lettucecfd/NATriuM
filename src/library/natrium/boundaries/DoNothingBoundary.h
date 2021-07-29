/*
 * SLFirstOrderBounceBack.h
 *
 *  Created on: 29.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_DONOTHINGBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_DONOTHINGBOUNDARY_H_

#include "../utilities/BasicNames.h"
#include "Boundary.h"
#include "BoundaryTools.h"

namespace natrium {

/**
 * @short simple first-order bounce-back for no-slip boundaries
 */
template<size_t dim>
class DoNothingBoundary: public Boundary<dim> {
public:
	/**
	 * @short Constructor
	 */
    DoNothingBoundary(size_t boundary_id) :
			Boundary<dim>(boundary_id, DO_NOTHING_BC, PrescribedBoundaryValues<dim>(dealii::Tensor<1,dim>())) {

	}

	/// destructor
	virtual ~DoNothingBoundary() {

	}

	/////////////////////////////////////////////////
	// FLAGS ////////////////////////////////////////
	/////////////////////////////////////////////////

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const{
		return false;
	}

	/** @short is the boundary a linear flux boundary as in SEDG-LBM (affine linear in the distributions f)?
	 */
	virtual bool isDGSupported() const{
		return false;
	}

	/** @short is the boundary set up to work with semi-Lagrangian streaming
	 */
	virtual bool isSLSupported() const{
		return true;
	}



	/**
	 * @short a function returning the values that are to be calculated at the boundaries
	 */
	virtual BoundaryFlags getUpdateFlags() const {
		return only_distributions;
	}

	/**
	 * @short calculate boundary values. This function overrides the virtual function of Boundary<dim>.
	 * @param fe_boundary_values An instance of FEBoundaryValues, that stores the flow variables at t-dt (here only the distribution functions)
	 * @param q_point the local index of the boundary hit point, i.e. its position in the fe_boundary_values
	 * @param destination defines the degree of freedom and discrete direction that the calculated value is assigned to
	 * 			(usually defines a point close to the boundary)
	 * @param eps Not used here. Defines the point in time at which the boundary is hit.
	 * @param t global time. Only relevant for time-dependent boundary conditions.
	 */
	virtual void calculateBoundaryValues(
			FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
			const LagrangianPathDestination& destination, double eps,
			double t) {

	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLFIRSTORDERBOUNCEBACK_H_ */
