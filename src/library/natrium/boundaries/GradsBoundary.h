/*
 * GradsBoundaryCondition.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_

#include <array>

#include "deal.II/base/tensor.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/function.h"

#include "GradsFunction.h"
#include "DoFBoundary.h"

#include "../advection/AdvectionOperator.h"
#include "../stencils/Stencil.h"
#include "../solver/DistributionFunctions.h"
#include "../utilities/BasicNames.h"

namespace natrium {

enum PrescribedQuantity {
	PRESCRIBED_VELOCITY, PRESCRIBED_PRESSURE
};

template<size_t dim, PrescribedQuantity prescribed_quantity>
class GradsBoundary: public DoFBoundary<dim> {
private:
	boost::shared_ptr<dealii::Function<dim> > m_boundaryValues;
public:
	/** @short This constructor assigns the Boundary condition with arbitrary density and velocity
	 *         to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] boundaryDensity A dealii::Function<dim> that defines the prescribed density at the boundary.
	 *  @param[in] boundaryVelocity A dealii::Function<dim> that defines the prescribed velocity at the boundary.
	 */
	GradsBoundary(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundary_values);

	virtual ~GradsBoundary() {

	}

	/**
	 * @short Apply the boundary condition (calculate unkown distributions)
	 */
	virtual void apply(DistributionFunctions& f, const distributed_vector& rho,
			const vector<distributed_vector>& u,
			const AdvectionOperator<dim>& advection, double beta,
			const Stencil& stencil);
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_ */
