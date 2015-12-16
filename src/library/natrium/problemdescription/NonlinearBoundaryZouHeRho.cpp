/*
 * NonlinearBoundaryZouHeRho.cpp
 *
 *  Created on: 15.12.2015
 *      Author: akraem3m
 */

#include "NonlinearBoundaryZouHeRho.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include "../advection/AdvectionOperator.h"
#include "../stencils/Stencil.h"

namespace natrium {

template<size_t dim>
NonlinearBoundaryZouHeRho<dim>::NonlinearBoundaryZouHeRho(
		boost::shared_ptr<dealii::Function<dim> > boundary_pressure,
		size_t direction) :
		m_boundaryPressure(boundary_pressure), m_direction(direction) {
	m_outwardNormal = dealii::Tensor<1, dim>();
	m_outwardNormal(dealii::GeometryInfo<dim>::unit_normal_direction[direction]) =
			dealii::GeometryInfo<dim>::unit_normal_orientation[direction];
	m_sign = dealii::GeometryInfo<dim>::unit_normal_orientation[direction];
}

template<size_t dim>
void NonlinearBoundaryZouHeRho<dim>::updateNonlinearBoundaryValues() const {

	// get all required data
	distributed_block_vector& boundary_vector =
			NonlinearBoundary<dim>::getBoundaryVector();
	const DistributionFunctions& f = NonlinearBoundary<dim>::getF();
	const distributed_vector& rho = NonlinearBoundary<dim>::getRho();
	const vector<distributed_vector>& u = NonlinearBoundary<dim>::getU();
	const dealii::UpdateFlags& update_flags =
			NonlinearBoundary<dim>::getUpdateFlags();
	const dealii::FE_DGQArbitraryNodes<dim>& fe =
			*(NonlinearBoundary<dim>::getAdvectionOperator().getFe());
	const dealii::QGaussLobatto<dim> & qGLL =
			*(NonlinearBoundary<dim>::getAdvectionOperator().getQuadrature());
	const dealii::DoFHandler<dim>& dof_handler =
			*(NonlinearBoundary<dim>::getAdvectionOperator().getDoFHandler());
	const dealii::MappingQ<dim>& mapping =
			NonlinearBoundary<dim>::getAdvectionOperator().getMapping();
	dealii::FEValues<dim> fe_values(mapping, fe, qGLL, update_flags);
	size_t this_boundary_id = NonlinearBoundary<dim>::getBoundaryIndicator();
	const Stencil& stencil = NonlinearBoundary<dim>::getStencil();

	const size_t dofs_per_cell = fe.dofs_per_cell;
	const size_t faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
	std::vector<dealii::types::global_dof_index> localDoFIndices(dofs_per_cell);
	const vector<std::map<size_t, size_t> >& q_index_to_facedof =
			*(NonlinearBoundary<dim>::getAdvectionOperator().getQIndexToFacedof());

	const vector<numeric_vector>& e = stencil.getDirections();
	const double scaling = stencil.getScaling();
	const size_t Q = stencil.getQ();

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
			const std::vector<dealii::Point<dim> > & points =
					fe_values.get_quadrature_points();
			const vector<dealii::Tensor<1, dim> > &normals =
					fe_values.get_all_normal_vectors();

			// loop over all quadrature points at the face
			for (size_t q = 0; q < fe_values.n_quadrature_points; q++) {
				size_t this_dof = q_index_to_facedof.at(i).at(q);
				assert(fe_values.shape_value(this_dof, q) > 0);

				// calculate velocity
				double p_in = m_boundaryPressure(points.at(this_dof));
				double rho_in = 1 + p_in / stencil.getSpeedOfSoundSquare();
				double u = 0;
				for (size_t j = 0; j < Q; j++) {
					size_t exn =
							(e.at(j)(
									dealii::GeometryInfo<dim>::unit_normal_direction[m_direction])
									* m_sign);
					if (exn > 1e-15) {
						u += 2 * f.at(j)(this_dof);
					} else if (exn > -1e-15) {
						//n.ei == 0
						u += f.at(j)(this_dof);
					}
					u += 1.0;
					u *= (-scaling);
					u *= m_sign;
				}

				// calculate fluxes (f_alpha - f_alpha^+)
				// as a bounce-back of non-equilibrium-distributions
				double prefactor = JxW.at(q);
				for (size_t j = 0; j < Q; j++) {
					size_t exn =
							(e.at(j)(
									dealii::GeometryInfo<dim>::unit_normal_direction[m_direction])
									* m_sign);
					if (exn < 0) {
						size_t opposite = stencil.getIndexOfOppositeDirection(
								j);
						// calculate scalar product
						double exu =
								u
										* e.at(j)(
												dealii::GeometryInfo<dim>::unit_normal_direction[m_direction]);
						double flux = prefactor * exn
								* (f.at(j)(this_dof) - f.at(opposite)(this_dof))
								- 2 * prefactor * stencil.getWeight(j) * rho_in
										* exu / stencil.getSpeedOfSoundSquare();
						boundary_vector.block(j - 1)(this_dof) += flux;
					} /* if exn < 0 */
				} /* for all directions */
			} /* for local quadrature points / dofs */
		} /* for all faces */
	} /* for all cells */
}

} /* namespace natrium */
