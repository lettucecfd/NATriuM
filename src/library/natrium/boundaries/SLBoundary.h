/*
 * SLBoundary.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_

#include "Boundary.h"

#include "deal.II/base/function.h"
#include "deal.II/base/tensor_function.h"

#include "../solver/DistributionFunctions.h"
#include "../stencils/Stencil.h"
#include "../advection/SemiLagrangianTools.h"
#include "../utilities/BasicNames.h"
#include "../advection/AdvectionOperator.h"
#include "BoundaryFlags.h"
#include "BoundaryTools.h"
#include "FEBoundaryValues.h"

namespace natrium {


template<size_t dim>
class SLBoundary: public Boundary<dim> {
private:
	size_t m_boundaryIndicator;
	PrescribedQuantities<dim> m_prescribedQuantities;
public:

	SLBoundary(size_t boundaryIndicator,
			const PrescribedQuantities<dim>& quantities) :
			m_boundaryIndicator(boundaryIndicator), m_prescribedQuantities(
					quantities) {

	}

	virtual ~SLBoundary() {
	}
	;

	virtual bool isPeriodic() const {
		return false;
	}
	virtual bool isLinearFluxBoundary() const {
		return false;
	}
	virtual bool isSLBoundary() const {
		return true;
	}

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}

	/**
	 * @short calculate boundary values. This pure virtual function has to be overriden by derived classes.
	 * @param fe_boundary_values An instance of FEBoundaryValues, that stores the flow variables at t-dt (here only the distribution functions)
	 * @param q_point the local index of the boundary hit point, i.e. its position in the fe_boundary_values
	 * @param destination defines the degree of freedom and discrete direction that the calculated value is assigned to
	 * 			(usually defines a point close to the boundary)
	 * @param eps  Defines the point in time at which the boundary is hit.
	 * @param t global time. Only relevant for time-dependent boundary conditions.
	 */
	virtual void calculateBoundaryValues(FEBoundaryValues<dim>& fe_boundary_values,
				size_t q_point, const LagrangianPathDestination& destination,
				double eps, double t) = 0;

	/**
	 * @short Get update flags. This pure virtual function has to be overriden by derived classes.
	 * @return the flags specifying the required values that have to be calculated from the flow field at t-dt,
	 * 			e.g. only_distributions, boundary_rho, boundary_drho_dt, ...
	 */
	virtual BoundaryFlags getUpdateFlags() const = 0;

	const PrescribedQuantities<dim>& getPrescribedQuantities() const {
		return m_prescribedQuantities;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_ */

