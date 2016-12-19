/*
 * SLFirstOrderBounceBack.h
 *
 *  Created on: 29.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLFIRSTORDERBOUNCEBACK_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLFIRSTORDERBOUNCEBACK_H_

#include "../utilities/BasicNames.h"
#include "SLBoundary.h"
#include "BoundaryTools.h"

namespace natrium {

/**
 * @short simple first-order bounce-back for no-slip boundaries
 */
template<size_t dim>
class SLFirstOrderBounceBack: public SLBoundary<dim> {
public:
	/**
	 * @short Constructor
	 */
	SLFirstOrderBounceBack(size_t boundary_id) :
			SLBoundary<dim>(boundary_id, PrescribedQuantities<dim>()) {

	}

	/// destructor
	virtual ~SLFirstOrderBounceBack() {

	}

	/**
	 * @short a function returning the values that are to be calculated at the boundaries
	 */
	virtual BoundaryFlags getUpdateFlags() const {
		return only_distributions;
	}

	/**
	 * @short calculate boundary values. This function overrides the virtual function of SLBoundary.
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
		// get opposite direction
		size_t opposite_i =
				fe_boundary_values.getData().m_stencil.getIndexOfOppositeDirection(
						destination.direction);
		// get opposite distribution
		double opposite_fi = fe_boundary_values.getDistribution(opposite_i,
				q_point);
		// assign opposite distribution to destination dof

		// does only work for points with support on boundary, thus:
		if (eps > 1e-7) {
			LOG(WARNING) << "Epsilon was " << eps
					<< " in SLFirstOrderBounceBack::calculateBoundaryValues(),"
							"which indicates that the destination point of the Lagrangian path ("
					<< fe_boundary_values.getPoint(q_point)
					<< ") was not at the boundary. That might lead to unwanted behavior."
					<< "If you really want to use the 1st order bounce-back scheme, you should consider "
							"having smaller time steps (CFL < 1 recommended)."
					<< endl;
			//throw BoundaryException()
		}
		//fe_boundary_values.getData().m_fnew.at(destination.direction)(destination.index) = opposite_fi;
		fe_boundary_values.getData().m_fnew.at(destination.direction)(
				destination.index) = fe_boundary_values.getData().m_fnew.at(
				opposite_i)(destination.index);
	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLFIRSTORDERBOUNCEBACK_H_ */
