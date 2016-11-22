/*
 * GradsBoundaryCondition.cpp
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#include "GradsBoundary.h"

namespace natrium {

template<size_t dim>
void GradsBoundary<dim>::updateMacroscopic(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
		const dealii::FEValues<dim>& fe_values, size_t q_point) {

	// check that object has been initialized properly via reinit
	assert(g.m_viscosity > 0);
	assert(g.m_dt > 0); // time step

	// check function parameters
	assert(fe_values.get_fe().n_dofs_per_cell() == b.local_dofs.size());
	assert (b.dofs_per_cell == b.local_dofs.size());

	// initalize macroscopic variables and dof variables
	b.rho = 0;
	b.rho_u = 0;
	b.u = 0;
	for (size_t i = 0; i < b.dofs_per_cell; i++) {
		b.rho_dof.at(i) = 0;
		b.rho_u_dof.at(i) = 0;
		b.u_dof.at(i) = 0;
	}

	// calculate macroscopic variables and dof variables
	for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
		for (size_t qi = 0; qi < g.m_Q; qi++) {
			b.rho_dof.at(dof) += g.m_fold.at(qi)(b.local_dofs.at(dof));
			for (size_t i = 0; i < dim; i++) {
				b.rho_u_dof.at(dof)[i] += g.m_fold.at(qi)(b.local_dofs.at(dof))
						* g.m_stencil.getDirection(qi)[i];
			}
		}
		b.u_dof.at(dof) = b.rho_u_dof.at(dof);
		b.u_dof.at(dof) /= b.rho_dof.at(dof);
		b.rho += b.rho_dof.at(dof) * fe_values.shape_value(dof, q_point);
		b.rho_u += b.rho_u_dof.at(dof) * fe_values.shape_value(dof, q_point);
		b.u += b.u_dof.at(dof) * fe_values.shape_value(dof, q_point);
	}
	assert((b.u * b.rho - b.rho_u).norm() < 1e-10);
}

template<size_t dim>
void GradsBoundary<dim>::calculateWallValues(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
		const dealii::FEValues<dim>& fe_values, size_t q_point) {
	b.rho_w = 0;
	b.drhodt_w = 0;
	b.u_w = 0;
	b.dudt_w = 0;
	b.dudx_w = 0;
	b.d2udxdt_w = 0;

	// calculate rho_w: first divergence, then rest
	for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
		for (size_t i = 0; i < dim; i++) {
			// derivative drhoui_dxi
			b.drhodt_w += fe_values.shape_grad(dof, q_point)[i]
					* b.rho_u_dof.at(dof)[i];
		}
	}
	b.drhodt_w *= -1;

	// calculate u_w: first viscous term, then derivative part, then rest
	for (size_t i = 0; i < dim; i++) {
		// viscous term
		for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
			for (size_t k = 0; k < dim; k++) {
				b.dudt_w[i] += fe_values.shape_hessian(dof, q_point)[k][k]
						* b.u_dof.at(dof)[i];
				b.dudt_w[i] += fe_values.shape_hessian(dof, q_point)[i][k]
						* b.u_dof.at(dof)[k];
			}
		}
		b.dudt_w[i] *= (-g.m_viscosity);
		for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
			// nonlinear term
			for (size_t k = 0; k < dim; k++) {
				b.dudt_w[i] += b.u[k] * fe_values.shape_grad(dof, q_point)[k]
						* b.u_dof.at(dof)[i];
			}
			// pressure term
			b.dudt_w[i] += g.m_cs2 / b.rho * fe_values.shape_grad(dof, q_point)[i]
					* b.rho_dof.at(dof);
		}
		b.u_w[i] *= -1;
	}

	// calculate dudx_w
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			// viscous term
			for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
				for (size_t k = 0; k < dim; k++) {
					b.d2udxdt_w[i][j] += fe_values.shape_3rd_derivative(dof,
							q_point)[k][k][j] * b.u_dof.at(dof)[i];
					b.d2udxdt_w[i][j] += fe_values.shape_3rd_derivative(dof,
							q_point)[i][k][j] * b.u_dof.at(dof)[k];
				}
			}
			b.d2udxdt_w[i][j] *= (-g.m_viscosity);
			for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
				for (size_t k = 0; k < dim; k++) {
					// nonlinear term 1
					b.d2udxdt_w[i][j] += fe_values.shape_grad(dof, q_point)[j]
							* b.u_dof.at(dof)[k]
							* fe_values.shape_grad(dof, q_point)[k]
							* b.u_dof.at(dof)[i];
					// nonlinear term 2
					b.d2udxdt_w[i][j] += b.u[j]
							* fe_values.shape_hessian(dof, q_point)[j][k]
							* b.u_dof.at(dof)[i];

				}
				// pressure term 1
				b.d2udxdt_w[i][j] += g.m_cs2 / (b.rho * b.rho)
						* fe_values.shape_hessian(dof, q_point)[i][j]
						* b.rho_dof.at(dof);
				// pressure term 2
				b.d2udxdt_w[i][j] += g.m_cs2 / (b.rho * b.rho)
						* fe_values.shape_grad(dof, q_point)[i]
						* b.rho_dof.at(dof)
						* fe_values.shape_grad(dof, q_point)[j]
						* b.rho_dof.at(dof);
			}
			b.d2udxdt_w[i] *= -1;
			// + dudx_w, which is not stored, but calculated by:
			for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
				b.dudx_w[i][j] += fe_values.shape_grad(dof, q_point)[j]
						* b.u_dof.at(dof)[i];
			}
		}
	}
}
template<size_t dim>
void GradsBoundary<dim>::applyWallValues(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
		const dealii::FEValues<dim>& fe_values, size_t q_point,
		const LagrangianPathDestination& destination, double eps) {

	assert(eps > 0); // local time to boundary
	assert(g.m_Q > 0);


	// calculate Grad's approximation of the distribution functions
	Tensor < 2, dim > P;
	Tensor < 2, dim > dudx_w_eps = b.dudx_w + (g.m_dt  - eps) * b.d2udxdt_w;
	Tensor < 1, dim > u_w_eps = b.u_w +  (g.m_dt  - eps) * b.dudt_w;
	double rho_w_eps =  b.rho_w + (g.m_dt  - eps) * b.drhodt_w;
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			P[i][j] = rho_w_eps * g.m_cs2 * (i == j) + rho_w_eps * u_w_eps[i] * u_w_eps[j]
					- rho_w_eps * (g.m_viscosity + g.m_cs2 / 2.0)
							* (dudx_w_eps[i][j] + dudx_w_eps[j][i]);
		}
	}

	vector<double> f;
	f.resize(g.m_Q);
	GradsFunction<dim>(f, g.m_stencil, b.rho_w + (g.m_dt - eps) * b.drhodt_w, (b.u_w + (g.m_dt -eps) * b.dudt_w) * b.rho_w, P);
	// TODO: calculate here

	g.m_fnew.at(destination.direction)(destination.index) = f.at(
			destination.direction);


}

template<size_t dim>
void GradsBoundary<dim>::calculateBoundaryValues(const GlobalBoundaryData& g, LocalBoundaryData<dim>& b,
		const dealii::FEValues<dim>& fe_values, size_t q_point,
		const LagrangianPathDestination& destination, double eps) {

	// TODO change so that this gets an instance of BoundaryHitsAtCell
	updateMacroscopic(g,b, fe_values, q_point);

	calculateWallValues(g,b, fe_values, q_point);

	applyWallValues(g,b, fe_values, q_point, destination, eps);

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

template class GradsBoundary<2> ;
template class GradsBoundary<3> ;

}
/* namespace natrium */
