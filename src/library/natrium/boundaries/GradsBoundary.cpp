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
		boost::shared_ptr<dealii::Function<dim> > boundary_values):
DoFBoundary<dim>(boundaryIndicator),
m_boundaryValues(boundary_values) {

}

template<size_t dim, PrescribedQuantity prescribed_quantity>
void GradsBoundary<dim, prescribed_quantity>::apply(DistributionFunctions& f,
		const distributed_block_vector rho,
		const vector<distributed_block_vector>& u,
		const AdvectionOperator<dim>& advection, double beta,
		const Stencil& stencil) {

	const dealii::DoFHandler<dim>& dof_handler = *advection.getDoFHandler();
	const size_t dofs_per_cell = advection.getNumberOfDoFsPerCell();
	const size_t faces_per_cell = dealii::GeometryInfo < dim > ::faces_per_cell;
	std::vector < dealii::types::global_dof_index
	> local_dof_indices(dofs_per_cell);
	const vector<std::map<size_t, size_t> >& q_index_to_facedof = advection.getQIndexToFacedof();
	const dealii::UpdateFlags updateFlags = dealii::update_gradients
	| dealii::update_quadrature_points | dealii::update_normal_vectors;

	dealii::FEFaceValues < dim
	> fe_face_values(advection.getMapping(), *advection.getFe(),
			*advection.getFaceQuadrature(), updateFlags);

	/*std::vector<dealii::types::global_dof_index> neighbor_indices(
	 dofs_per_cell);
	 dealii::FEValues<dim> neighbor_fe_values(advection.getMapping(),
	 *advection.getFe(), *advection.getQuadrature(), updateFlags);
	 */

	double rho_tgt;
	//double rho_bb;
	//double rho_s;
	dealii::Tensor < 1, dim > j_tgt;
	dealii::Tensor < 2, dim > P;
	dealii::Tensor < 2, dim > P_eq;
	dealii::Tensor < 2, dim > P_neq;
	vector<double> f_grad;
	f_grad.resize(stencil.getQ());
	std::vector < std::vector<Tensor<1, dim> > > u_gradients; // .at(i).at(j)[k] denotes du_i / dx_k (x_j)

	// loop over all cells
	dealii::Tensor < 1, dim > density_gradient;
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

			cell->face(fc)->get_dof_indices(local_dof_indices);

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
					j_tgt[j] = this->getBoundaryVelocity()->value(
							fe_face_values.quadrature_point(q), j);

					for (size_t k = 0; k < dim; k++) {
						P_eq[j][k] += u.at(j)(local_dof_indices.at(i))
						* u.at(j)(local_dof_indices.at(i));
						P_neq[j][k] = u_gradients.at(j).at(i)[k]
						+ u_gradients.at(k).at(i)[j];
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
					if (vectorToTensor < dim
							> (stencil.getDirection(alpha))
							* fe_face_values.normal_vector(i) < 0) {
						f.at(alpha)(local_dof_indices.at(i)) = f_grad.at(alpha);
					}
				}

			}
		}

	}
}

template class GradsBoundary<2, PRESCRIBED_VELOCITY> ;
template class GradsBoundary<3, PRESCRIBED_VELOCITY> ;

}
/* namespace natrium */
