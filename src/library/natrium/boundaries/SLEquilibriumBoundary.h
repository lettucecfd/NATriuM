/*
 * SLEquilibriumBoundary.h
 *
 *  Created on: 28.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_

#include "../utilities/BasicNames.h"
#include "SLBoundary.h"

namespace natrium {

/**
 * @short Boundary condition that sets an equilibrium with respect to a prescribed velocity or pressure.
 */
template<size_t dim>
class SLEquilibriumBoundary: public SLBoundary<dim> {
public:
	SLEquilibriumBoundary(size_t boundary_id, dealii::Tensor<1, dim>& velocity) :
			SLBoundary<dim>(boundary_id, PrescribedQuantities<dim>(velocity)) {
	}

	virtual ~SLEquilibriumBoundary();

	virtual BoundaryFlags getUpdateFlags() const {
		BoundaryFlags flags = only_distributions;
		if (boundary_u & SLBoundary<dim>::getPrescribedQuantities()) {
			flags |= (boundary_rho | boundary_drho_dt);

		}
		if ((boundary_rho & SLBoundary<dim>::getPrescribedQuantities())
				or (boundary_p & SLBoundary<dim>::getPrescribedQuantities())) {
			flags |= (boundary_u | boundary_du_dt);
		}
		return flags;
	}

	virtual void calculateBoundaryValues(
			FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
			const LagrangianPathDestination& destination, double eps,
			double t) {

		// get density
		double rho;
		if (boundary_p & SLBoundary<dim>::getPrescribedQuantities()) {
			// set time of this function object for time-dependent boundary conditions
			SLBoundary<dim>::getPrescribedQuantities().getPressure().set_time(
					t - eps);
			// evaluate boundary conditions
			double p =
					SLBoundary<dim>::getPrescribedQuantities().getPressure().value(
							fe_boundary_values.getPoint(q_point));
			// obtain density by ideal gas law
			rho = p / fe_boundary_values.getData().m_cs2;
		} else {
			// calculate density at t - eps by Taylor expansion
			rho = fe_boundary_values.getRho()
					+ (fe_boundary_values.getData().m_dt - eps)
							* fe_boundary_values.getDrhodt();
		}

		// get velocity
		dealii::Tensor<1, dim> u;
		if (boundary_u & m_prescribedQuantities) {
			// set time of this function object for time-dependent boundary conditions
			SLBoundary<dim>::getPrescribedQuantities().getVelocity().set_time(
					t - eps);
			// evaluate boundary conditions
			dealii::Vector<double> tmp(dim);
			SLBoundary<dim>::getPrescribedQuantities().getVelocity().vector_value(
					tmp, fe_boundary_values.getPoint(q_point));
			// vector to tensor
			for (size_t i = 0; i < dim; i++) {
				u[i] = tmp(i);
			}
		} else {
			// Taylor expansion
			u = fe_boundary_values.getU()
					+ (fe_boundary_values.getData().m_dt - eps)
							* fe_boundary_values.getDudt();
		}

		// calculate BGK equilibrium
		dealii::Tensor<1, dim> e_i = vectorToTensor<dim>(
				fe_boundary_values.getData().m_stencil.getWeight(
						destination.direction));
		double w_i = fe_boundary_values.getData().m_stencil.getWeight(
				destination.direction);
		double cs2 = fe_boundary_values.getData().m_cs2;
		fe_boundary_values.getData().m_fnew.at(destination.direction)(
				destination.index) = rho * w_i
				* (1 + (e_i * u) / cs2 + pow(e_i * u, 2) / (2 * cs2 * cs2)
						- (u * u) / (2 * cs2));

	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_ */
