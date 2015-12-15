/*
 * NonlinearBoundaryZouHeRho.cpp
 *
 *  Created on: 15.12.2015
 *      Author: akraem3m
 */

#include "NonlinearBoundaryZouHeRho.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"

namespace natrium {

template<size_t dim>
NonlinearBoundaryZouHeRho<dim>::NonlinearBoundaryZouHeRho(
		boost::shared_ptr<dealii::Function<dim> > boundary_density) :
		m_boundaryDensity(boundary_density) {
}

template<size_t dim>
void NonlinearBoundaryZouHeRho<dim>::updateNonlinearBoundaryValues() const {

	// get all required data
	distributed_vector& boundary_vector = getBoundaryVector();
	const DistributionFunctions& f = getF();
	const distributed_vector& rho = getRho();
	const vector<distributed_vector>& u = getU();
	const dealii::UpdateFlags& update_flags = getUpdateFlags();
	const dealii::FE_DGQArbitraryNodes<dim>& fe =
			*(getAdvectionOperator().getFe());
	const dealii::QGaussLobatto<dim> & qGLL =
			*(getAdvectionOperator().getQuadrature());
	const dealii::DoFHandler<dim>& dof_handler =
			*(getAdvectionOperator().getDoFHandler());
	const dealii::MappingQ<dim>& mapping = getAdvectionOperator().getMapping();
	dealii::FEValues<dim> fe_values(mapping, fe, qGLL, update_flags);
	size_t this_boundary_id = getBoundaryIndicator();

	const size_t dofs_per_cell = fe.dofs_per_cell;
	const size_t faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
	std::vector<dealii::types::global_dof_index> localDoFIndices(dofs_per_cell);
	const vector<std::map<size_t, size_t> >& q_index_to_facedof = *(getAdvectionOperator().getQIndexToFacedof());

	//////////////////////////////////
	// LOOP OVER ALL BOUNDARY FACES //
	//////////////////////////////////
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler->begin_active(), endc = dof_handler->end();
	for (; cell != endc; ++cell) {
		if (!cell->is_locally_owned()) {
			continue;
		}
		for (size_t i = 0; i < faces_per_cell; i++) {
			//Faces at boundary
			if (!cell->face(i)->at_boundary()) {
				continue;
			}
			size_t boundary_id = cell->face(i)->boundary_id();
			if (boundary_id != this_boundary_id) {
				continue;
			}
			// calculate the fe values for the cell
			fe_values.reinit(cell, i);
			cell->get_dof_indices(localDoFIndices);

			// get data
			const vector<double> &JxW = fe_values.get_JxW_values();
			const vector<dealii::Tensor<1, dim> > &normals =
					fe_values.get_all_normal_vectors();
			const dealii::Vector& e_alpha = getStencil()->getDirection(1);

			// loop over all quadrature points at the face
			for (size_t q = 0; q < fe_values.n_quadrature_points; q++) {
				size_t thisDoF = q_index_to_facedof.at(i).at(q);
				assert(fe_values.shape_value(thisDoF, q) > 0);
			}

			// calculate fluxes

		} /* for all faces */
	}/* for all cells */

}

} /* namespace natrium */
