/*
 * GradsBoundaryCondition.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_

//#include "SLBoundary.h"
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

//template<size_t dim>
//class GradsBoundary: public SLBoundary<dim> {
//public:
//
//	GradsBoundary(size_t boundary_id,
//			const PrescribedBoundaryValues<dim>& quantities) :
//			SLBoundary<dim>(boundary_id, quantities) {
//
//	}
//
//	GradsBoundary(size_t boundary_id, double pressure) :
//			SLBoundary<dim>(boundary_id, PrescribedBoundaryValues<dim>(pressure)) {
//
//	}
//
//	GradsBoundary(size_t boundary_id,
//			boost::shared_ptr<dealii::TensorFunction<2, dim> > function) :
//			SLBoundary<dim>(boundary_id, PrescribedBoundaryValues<dim>(function)) {
//
//	}
//
//	GradsBoundary(size_t boundary_id, dealii::Tensor<1, dim>& velocity) :
//			SLBoundary<dim>(boundary_id, PrescribedBoundaryValues<dim>(velocity)) {
//
//	}
//
//	GradsBoundary(size_t boundary_id, dealii::Tensor<2, dim>& velocity_gradient) :
//			SLBoundary<dim>(boundary_id,
//					PrescribedBoundaryValues<dim>(velocity_gradient)) {
//
//	}
//
//	GradsBoundary(size_t boundary_id, const dealii::Vector<double>& velocity) :
//			SLBoundary<dim>(boundary_id, PrescribedBoundaryValues<dim>(velocity)) {
//	}
//
//	virtual ~GradsBoundary() {
//	}
//
//	virtual BoundaryFlags getUpdateFlags() const {
//		return only_distributions;
//	}
//
//	virtual void calculateBoundaryValues(FEBoundaryValues<dim>& fe_boundary_values,
//				size_t q_point, const LagrangianPathDestination& destination,
//				double eps, double t){
//
//	}
//
//
//};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_ */
