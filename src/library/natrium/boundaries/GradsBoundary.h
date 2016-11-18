/*
 * GradsBoundaryCondition.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_

#include "SLBoundary.h"
#include <array>

#include "deal.II/base/tensor.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/function.h"

#include "GradsFunction.h"
#include "../advection/AdvectionOperator.h"
#include "../stencils/Stencil.h"
#include "../solver/DistributionFunctions.h"
#include "../utilities/BasicNames.h"

namespace natrium {

enum PrescribedQuantity {
	PRESCRIBED_VELOCITY, PRESCRIBED_PRESSURE
};

template<size_t dim, PrescribedQuantity prescribed_quantity>
class GradsBoundary: public SLBoundary<dim> {
private:
	boost::shared_ptr<Stencil> m_stencil;
	double m_viscosity;
	double m_dt;
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

	GradsBoundary(size_t boundaryIndicator, double pressure);

	//GradsBoundary(size_t boundaryIndicator,
	//		dealii::Tensor<1,dim> velocity);
	/// constructor
	GradsBoundary(size_t boundaryIndicator,
			const dealii::Vector<double>& velocity);

	virtual ~GradsBoundary() {

	}

	/**
	 * @short Apply the boundary condition (calculate unkown distributions)
	 */
	virtual void calculateBoundaryValues(const DistributionFunctions& f_old,
			DistributionFunctions& f_new,
			const std::vector<dealii::types::global_dof_index> & local_dofs,
			const dealii::FEValues<dim>&, size_t q_point,
			const LagrangianPathDestination& destination, double dt) const;

	virtual dealii::UpdateFlags getUpdateFlags() const {
		return dealii::update_values | dealii::update_gradients
				| dealii::update_3rd_derivatives | dealii::update_hessians;

	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_ */
