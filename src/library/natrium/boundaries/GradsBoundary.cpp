/*
 * GradsBoundaryCondition.cpp
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#include "GradsBoundary.h"

namespace natrium {

template<size_t dim, PrescribedQuantity prescribed_quantity>
GradsBoundary<dim, prescribed_quantity>::GradsBoundary(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundary_values) :
		SLBoundary<dim>(boundaryIndicator), m_boundaryValues(boundary_values) {
	// TODO: assert n_components = dim

}

template<size_t dim, PrescribedQuantity prescribed_quantity>
GradsBoundary<dim, prescribed_quantity>::GradsBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity) :
		SLBoundary<dim>(boundaryIndicator), m_boundaryValues(
				boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(
						velocity)) {
	assert(prescribed_quantity == PRESCRIBED_VELOCITY);
}

template<size_t dim, PrescribedQuantity prescribed_quantity>
GradsBoundary<dim, prescribed_quantity>::GradsBoundary(size_t boundaryIndicator,
		double pressure) :
		SLBoundary<dim>(boundaryIndicator), m_boundaryValues(
				boost::make_shared<BoundaryTools::BoundaryPressure<dim> >(
						pressure)) {

	assert(prescribed_quantity == PRESCRIBED_PRESSURE);
}

template<size_t dim, PrescribedQuantity prescribed_quantity>
void GradsBoundary<dim, prescribed_quantity>::calculateBoundaryValues(
		const DistributionFunctions& f_old, DistributionFunctions& f_new,
		const std::vector<dealii::types::global_dof_index> & local_dofs,
		const dealii::FEValues<dim>& fe_values, size_t q_point,
		const LagrangianPathDestination& destination, double eps) const {

	// check that object has been initialized properly
	assert(m_viscosity > 0);
	assert(m_stencil);
	assert(m_dt > 0); // time step
	// check function parameters
	assert(eps > 0); // local time to boundary
	assert(fe_values.get_fe().n_dofs_per_cell() == local_dofs.size());

	// obtain constants
	const size_t dofs_per_cell = local_dofs.size();
	const size_t Q = m_stencil->getQ();
	const double cs2 = m_stencil->getSpeedOfSoundSquare();

	// initalize macroscopic variables and dof variables
	double rho = 0.0;
	Tensor < 1, dim > rho_u;
	Tensor < 1, dim > u;
	vector<double> rho_dof;
	rho_dof.resize(dofs_per_cell);
	vector<Tensor<1, dim> > rho_u_dof;
	rho_u_dof.resize(dofs_per_cell);
	vector<Tensor<1, dim> > u_dof;
	u_dof.resize(dofs_per_cell);

	// calculate macroscopic variables and dof variables
	for (size_t dof = 0; dof < dofs_per_cell; dof++) {
		for (size_t qi = 0; qi < Q; qi++) {
			rho_dof.at(dof) += f_old.at(qi)(local_dofs.at(dof));
			for (size_t i = 0; i < dim; i++) {
				rho_u_dof.at(dof)[i] += f_old.at(qi)(local_dofs.at(dof))
						* m_stencil->getDirection(qi)[i];
			}
		}
		u_dof.at(dof) = rho_u_dof.at(dof);
		u_dof.at(dof) /= rho_dof.at(dof);
		rho += rho_dof.at(dof) * fe_values.shape_value(dof, q_point);
		rho_u += rho_u_dof.at(dof) * fe_values.shape_value(dof, q_point);
		u += u_dof.at(dof) * fe_values.shape_value(dof, q_point);
	}
	assert((u * rho - rho_u).norm() < 1e-10);

	// calculate boundary values
	// TODO: incorporate prescribed values
	double rho_w = 0;
	Tensor < 1, dim > u_w;
	Tensor < 2, dim > dudx_w;

	// calculate rho_w: first divergence, then rest
	for (size_t dof = 0; dof < dofs_per_cell; dof++) {
		for (size_t i = 0; i < dim; i++) {
			// derivative drhoui_dxi
			rho_w += fe_values.shape_grad(dof, q_point)[i]
					* rho_u_dof.at(dof)[i];
		}
	}
	rho_w *= -(m_dt - eps);
	rho_w += rho;

	// calculate u_w: first viscous term, then derivative part, then rest
	for (size_t i = 0; i < dim; i++) {
		// viscous term
		for (size_t dof = 0; dof < dofs_per_cell; dof++) {
			for (size_t k = 0; k < dim; k++) {
				u_w[i] += fe_values.shape_hessian(dof, q_point)[k][k]
						* u_dof.at(dof)[i];
				u_w[i] += fe_values.shape_hessian(dof, q_point)[i][k]
						* u_dof.at(dof)[k];
			}
		}
		u_w[i] *= (-m_viscosity);
		for (size_t dof = 0; dof < dofs_per_cell; dof++) {
			// nonlinear term
			for (size_t k = 0; k < dim; k++) {
				u_w[i] += u[k] * fe_values.shape_grad(dof, q_point)[k]
						* u_dof.at(dof)[i];
			}
			// pressure term
			u_w[i] += cs2 / rho * fe_values.shape_grad(dof, q_point)[i]
					* rho_dof.at(dof);
		}
		u_w[i] *= -(m_dt - eps);
		u_w[i] += u[i];
	}

	// calculate dudx_w
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			// viscous term
			for (size_t dof = 0; dof < dofs_per_cell; dof++) {
				for (size_t k = 0; k < dim; k++) {
					dudx_w[i][j] +=
							fe_values.shape_3rd_derivative(dof, q_point)[k][k][j]
									* u_dof.at(dof)[i];
					dudx_w[i][j] +=
							fe_values.shape_3rd_derivative(dof, q_point)[i][k][j]
									* u_dof.at(dof)[k];
				}
			}
			dudx_w[i][j] *= (-m_viscosity);
			for (size_t dof = 0; dof < dofs_per_cell; dof++) {
				for (size_t k = 0; k < dim; k++) {
					// nonlinear term 1
					dudx_w[i][j] += fe_values.shape_grad(dof, q_point)[j]
							* u_dof.at(dof)[k]
							* fe_values.shape_grad(dof, q_point)[k]
							* u_dof.at(dof)[i];
					// nonlinear term 2
					dudx_w[i][j] += u[j]
							* fe_values.shape_hessian(dof, q_point)[j][k]
							* u_dof.at(dof)[i];

				}
				// pressure term 1
				dudx_w[i][j] += cs2 / (rho * rho)
						* fe_values.shape_hessian(dof, q_point)[i][j]
						* rho_dof.at(dof);
				// pressure term 2
				dudx_w[i][j] += cs2 / (rho * rho)
						* fe_values.shape_grad(dof, q_point)[i]
						* rho_dof.at(dof)
						* fe_values.shape_grad(dof, q_point)[j]
						* rho_dof.at(dof);
			}
			dudx_w[i] *= -(m_dt - eps);
			// + dudx_w, which is not stored, but calculated by:
			for (size_t dof = 0; dof < dofs_per_cell; dof++) {
				dudx_w[i][j] += fe_values.shape_grad(dof, q_point)[j]
						* u_dof.at(dof)[i];
			}
		}
	}



	/*
	 assert(u.size() == dim);
	 const dealii::DoFHandler<dim>& dof_handler = *advection.getDoFHandler();
	 const size_t dofs_per_cell = advection.getNumberOfDoFsPerCell();
	 const size_t faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
	 std::vector<dealii::types::global_dof_index> local_dof_indices(
	 dofs_per_cell);
	 const vector<std::map<size_t, size_t> >& q_index_to_facedof =
	 advection.getQIndexToFacedof();
	 assert(q_index_to_facedof.size() > 0);
	 const dealii::UpdateFlags updateFlags = dealii::update_gradients
	 | dealii::update_quadrature_points | dealii::update_normal_vectors;

	 dealii::FEFaceValues<dim> fe_face_values(advection.getMapping(),
	 *advection.getFe(), *advection.getFaceQuadrature(), updateFlags);

	 double rho_tgt;
	 //double rho_bb;
	 //double rho_s;
	 dealii::Tensor<1, dim> j_tgt;
	 dealii::Tensor<2, dim> P;
	 dealii::Tensor<2, dim> P_eq;
	 dealii::Tensor<2, dim> P_neq;
	 vector<double> f_grad;
	 f_grad.resize(stencil.getQ());
	 std::vector<std::vector<Tensor<1, dim> > > u_gradients(dim); // .at(i).at(j)[k] denotes du_i / dx_k (x_j)
	 for (size_t j = 0; j < dim; j++) {
	 u_gradients.at(j).resize(fe_face_values.n_quadrature_points);
	 }
	 // loop over all cells
	 dealii::Tensor<1, dim> density_gradient;
	 typename dealii::DoFHandler<dim>::active_cell_iterator cell =
	 dof_handler.begin_active(), endc = dof_handler.end();
	 for (; cell != endc; ++cell) {

	 if (!cell->is_locally_owned())
	 continue;

	 if (!cell->at_boundary())
	 continue;

	 // get global degrees of freedom

	 for (size_t fc = 0; fc < faces_per_cell; fc++) {

	 if (!cell->at_boundary(fc))
	 continue;
	 if (this->getBoundaryIndicator() != cell->face(fc)->boundary_id())
	 continue;
	 cell->get_dof_indices(local_dof_indices);

	 // calculate gradients of shape functions
	 fe_face_values.reinit(cell, fc);
	 for (size_t j = 0; j < dim; j++) {
	 fe_face_values.get_function_gradients(u.at(j),
	 u_gradients.at(j));
	 }

	 // for all dofs
	 for (size_t q = 0; q < fe_face_values.n_quadrature_points; q++) {
	 size_t i = q_index_to_facedof.at(fc).at(q);
	 P_eq = 0;
	 P_neq = 0;
	 j_tgt = 0;
	 rho_tgt = 1;
	 for (size_t j = 0; j < dim; j++) {
	 P_eq[j][j] = stencil.getSpeedOfSoundSquare();
	 j_tgt[j] = m_boundaryValues->value(
	 fe_face_values.quadrature_point(q), j);
	 for (size_t k = 0; k < dim; k++) {
	 ;
	 assert(local_dof_indices.size() > 0);
	 P_eq[j][k] += u.at(j)(local_dof_indices.at(i))
	 * u.at(j)(local_dof_indices.at(i));
	 P_neq[j][k] = u_gradients.at(j).at(q)[k]
	 + u_gradients.at(k).at(q)[j];
	 }
	 }
	 P_eq *= rho(local_dof_indices.at(i));
	 P_neq *= (-rho(local_dof_indices.at(i))
	 * stencil.getSpeedOfSoundSquare() / 2.0 / beta);
	 P = 0;
	 P += P_eq;
	 P += P_neq;

	 // TODO incorporate moving surfaces by using rho_tgt != 1
	 // TODO incorporate pressure boundaries by using rho_tgt != 1 and u_tgt != uw
	 GradsFunction<dim>(f_grad, stencil, rho_tgt, j_tgt, P);
	 for (size_t alpha = 0; alpha < stencil.getQ(); alpha++) {
	 if (vectorToTensor<dim>(stencil.getDirection(alpha))
	 * fe_face_values.normal_vector(q) < 0) {
	 f.at(alpha)(local_dof_indices.at(i)) = f_grad.at(alpha);
	 }
	 }

	 }
	 }

	 }
	 */
}

template class GradsBoundary<2, PRESCRIBED_VELOCITY> ;
template class GradsBoundary<3, PRESCRIBED_VELOCITY> ;

}
/* namespace natrium */
