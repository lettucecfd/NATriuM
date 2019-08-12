/*
 * SLEquilibriumBoundary.h
 *
 *  Created on: 28.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_

#include "../utilities/BasicNames.h"
#include "Boundary.h"
#include "BoundaryFlags.h"

namespace natrium {

/**
 * @short Boundary condition that sets an equilibrium with respect to a prescribed velocity or pressure.
 */
template<size_t dim>
class SLEquilibriumBoundary: public Boundary<dim> {
public:
	SLEquilibriumBoundary(size_t boundary_id, dealii::Tensor<1, dim>& velocity) :
			Boundary<dim>(boundary_id, VELOCITY_EQUILIBRIUM_BOUNDARY,
					PrescribedBoundaryValues<dim>(velocity)) {
	}

	//virtual ~SLEquilibriumBoundary();

	/////////////////////////////////////////////////
	// FLAGS ////////////////////////////////////////
	/////////////////////////////////////////////////

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const {
		return false;
	}

	/** @short is the boundary a linear flux boundary as in SEDG-LBM (affine linear in the distributions f)?
	 */
	virtual bool isDGSupported() const {
		return false;
	}

	/** @short is the boundary set up to work with semi-Lagrangian streaming
	 */
	virtual bool isSLSupported() const {
		return true;
	}

	virtual BoundaryFlags getUpdateFlags() const {
		BoundaryFlags flags = only_distributions;
		return flags;
	}

	virtual void calculateBoundaryValues(
			FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
			const LagrangianPathDestination& destination, double eps,
			double t) {

		const GlobalBoundaryData& data = fe_boundary_values.getData();
		const Stencil& stencil = data.m_stencil;

		// density unity
		double rho = 1.0;

		double feq = 0.0;
		double cs2 = stencil.getSpeedOfSoundSquare();

		// get velocity
		dealii::Tensor<1, dim> u;
		// set time of this function object for time-dependent boundary conditions

		assert(Boundary<dim>::getBoundaryValues().getVelocity());
		Boundary<dim>::getBoundaryValues().getVelocity()->set_time(t - eps);
		// evaluate boundary conditions
		dealii::Vector<double> tmp(dim);
		Boundary<dim>::getBoundaryValues().getVelocity()->vector_value(
				fe_boundary_values.getPoint(q_point), tmp);

		// vector to tensor
		for (size_t i = 0; i < dim; i++) {
			u[i] = tmp(i);
		}

		double eye[2][2] = { { 1, 0 }, { 0, 1 } };
		double uu_term = 0.0;
		for (size_t j = 0; j < 2; j++) {
			uu_term += -(u[j] * u[j]) / (2.0 * cs2);
		}

			double T0 = 1.0;
			double T1 = cs2 * (T0 - 1);

			double ue_term = 0.0;
			for (size_t j = 0; j < dim; j++) {
				ue_term += (u[j] * stencil.getDirection(destination.direction)[j]) / cs2;
			}
			feq = stencil.getWeight(destination.direction) * rho
					* (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
			for (int alp = 0; alp < 2; alp++) {
				for (int bet = 0; bet < 2; bet++) {
					feq += rho * stencil.getWeight(destination.direction) / (2.0 * cs2)
							* ((T0 - 1) * eye[alp][bet] * stencil.getDirection(destination.direction)[alp]
									* stencil.getDirection(destination.direction)[bet]
									- cs2 * eye[alp][bet]
											* (T0 - 1));
					for (int gam = 0; gam < 2; gam++) {

						feq +=
								stencil.getWeight(destination.direction) * rho
										/ (6. * cs2 * cs2 * cs2)
										* (u[alp] * u[bet]
												* u[gam]
												+ (T0 - 1) * cs2
														* (eye[alp][bet]
																* u[gam]
																+ eye[bet][gam]
																		* u[alp]
																+ eye[alp][gam]
																		* u[bet]))
										* (stencil.getDirection(destination.direction)[alp] * stencil.getDirection(destination.direction)[bet]
												* stencil.getDirection(destination.direction)[gam]
												- cs2
														* (stencil.getDirection(destination.direction)[gam]
																* eye[alp][bet]
																+ stencil.getDirection(destination.direction)[bet]
																		* eye[alp][gam]
																+ stencil.getDirection(destination.direction)[alp]
																		* eye[bet][gam]));

						for (int det = 0; det < 2; det++) {
							double power4 = stencil.getDirection(destination.direction)[alp] * stencil.getDirection(destination.direction)[bet]
									* stencil.getDirection(destination.direction)[gam] * stencil.getDirection(destination.direction)[det];
							double power2 = stencil.getDirection(destination.direction)[alp] * stencil.getDirection(destination.direction)[bet]
									* eye[gam][det]
									+ stencil.getDirection(destination.direction)[alp] * stencil.getDirection(destination.direction)[gam] * eye[bet][det]
									+ stencil.getDirection(destination.direction)[alp] * stencil.getDirection(destination.direction)[det] * eye[bet][gam]
									+ stencil.getDirection(destination.direction)[bet] * stencil.getDirection(destination.direction)[gam] * eye[alp][det]
									+ stencil.getDirection(destination.direction)[bet] * stencil.getDirection(destination.direction)[det] * eye[alp][gam]
									+ stencil.getDirection(destination.direction)[gam] * stencil.getDirection(destination.direction)[det] * eye[alp][bet];
							double power0 = eye[alp][bet] * eye[gam][det]
									+ eye[alp][gam] * eye[bet][det]
									+ eye[alp][det] * eye[bet][gam];
							double u4 = u[alp] * u[bet]
									* u[gam] * u[det];
							double u2 = u[alp] * u[bet]
									* eye[gam][det]
									+ u[alp] * u[gam]
											* eye[bet][det]
									+ u[alp] * u[det]
											* eye[bet][gam]
									+ u[bet] * u[gam]
											* eye[alp][det]
									+ u[bet] * u[det]
											* eye[alp][gam]
									+ u[gam] * u[det]
											* eye[alp][bet];
							double multieye = eye[alp][bet] * eye[gam][det]
									+ eye[alp][gam] * eye[bet][det]
									+ eye[alp][det] * eye[bet][gam];

							feq += stencil.getWeight(destination.direction) * rho
									/ (24. * cs2 * cs2 * cs2 * cs2)
									* (power4 - cs2 * power2
											+ cs2 * cs2 * power0)
									* (u4 + T1 * (u2 + T1 * multieye));

						}

					}
				}

			}


			fe_boundary_values.getData().m_fnew.at(destination.direction)(
								destination.index) = feq;

		/*
		 *
		 *
		 // get density
		 double rho;
		 if (Boundary<dim>::getBoundaryValues()) {
		 // set time of this function object for time-dependent boundary conditions
		 Boundary<dim>::getBoundaryValues().getPressure().set_time(
		 t - eps);
		 // evaluate boundary conditions
		 double p =
		 Boundary<dim>::getBoundaryValues().getPressure().value(
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
		 if (Boundary<dim>::getPrescribedQuantities()) {
		 // set time of this function object for time-dependent boundary conditions
		 Boundary<dim>::getBoundaryValues().getVelocity().set_time(
		 t - eps);
		 // evaluate boundary conditions
		 dealii::Vector<double> tmp(dim);
		 Boundary<dim>::getPrescribedValues().getVelocity().vector_value(
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
		 vector<numeric_vector> e_i(fe_boundary_values.getData().m_stencil.getDirections());

		 //dealii::Tensor<1, dim> e_i = vectorToTensor<dim>(
		 //		fe_boundary_values.getData().m_stencil.getWeight(
		 //				destination.direction));

		 double w_i = fe_boundary_values.getData().m_stencil.getWeight(
		 destination.direction);
		 double cs2 = fe_boundary_values.getData().m_cs2;
		 */

	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_ */
