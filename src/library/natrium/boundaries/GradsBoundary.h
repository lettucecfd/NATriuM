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

template<size_t dim>
class GradsBoundary: public SLBoundary<dim> {
public:

	GradsBoundary(size_t boundary_id,
			const PrescribedQuantities<dim>& quantities) :
			SLBoundary<dim>(boundary_id, quantities) {

	}

	GradsBoundary(size_t boundary_id, double pressure) :
			SLBoundary<dim>(boundary_id, PrescribedQuantities<dim>(pressure)) {

	}

	GradsBoundary(size_t boundary_id,
			boost::shared_ptr<dealii::TensorFunction<2, dim> > function) :
			SLBoundary<dim>(boundary_id, PrescribedQuantities<dim>(function)) {

	}

	GradsBoundary(size_t boundary_id, dealii::Tensor<1, dim>& velocity) :
			SLBoundary<dim>(boundary_id, PrescribedQuantities<dim>(velocity)) {

	}

	GradsBoundary(size_t boundary_id, dealii::Tensor<2, dim>& velocity_gradient) :
			SLBoundary<dim>(boundary_id,
					PrescribedQuantities<dim>(velocity_gradient)) {

	}

	GradsBoundary(size_t boundary_id, const dealii::Vector<double>& velocity) :
			SLBoundary<dim>(boundary_id, PrescribedQuantities<dim>(velocity)) {
	}

	virtual ~GradsBoundary() {
	}

	void updateMacroscopic(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
			const dealii::FEValues<dim>& fe_values, size_t q_point);

	void calculateWallValues(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
			const dealii::FEValues<dim>& fe_values, size_t q_point);

	void applyWallValues(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
			const dealii::FEValues<dim>& fe_values, size_t q_point,
			const LagrangianPathDestination& destination, double dt);

	/**
	 * @short Apply the boundary condition (calculate unkown distributions)
	 */
	virtual void calculateBoundaryValues(const GlobalBoundaryData& g,
			LocalBoundaryData<dim>& b, const dealii::FEValues<dim>& fe_values,
			size_t q_point, const LagrangianPathDestination& destination,
			double dt);

	virtual dealii::UpdateFlags getUpdateFlags() const {
		return dealii::update_values | dealii::update_gradients
				| dealii::update_3rd_derivatives | dealii::update_hessians;

	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_ */
