/*
 * DealiiExtensions.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 */

#include <utilities/DealiiExtensions.h>

#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/multigrid/mg_dof_handler.h>

#include <algorithm>
#include <numeric>

using namespace dealii;
using namespace DoFTools;

namespace natrium {

namespace DealIIExtensions {

template<class DH, class SparsityPattern>
void make_sparser_flux_sparsity_pattern(const DH &dof,
		SparsityPattern &sparsity, const ConstraintMatrix &constraints,
		const bool keep_constrained_dofs,
		const types::subdomain_id subdomain_id)

		{

	const types::global_dof_index n_dofs = dof.n_dofs();

	AssertDimension(sparsity.n_rows(), n_dofs);
	AssertDimension(sparsity.n_cols(), n_dofs);

	// If we have a distributed::Triangulation only allow locally_owned
	// subdomain. Not setting a subdomain is also okay, because we skip
	// ghost cells in the loop below.
	Assert(
			(dof.get_tria().locally_owned_subdomain() == numbers::invalid_subdomain_id) || (subdomain_id == numbers::invalid_subdomain_id) || (subdomain_id == dof.get_tria().locally_owned_subdomain()),
			ExcMessage ("For parallel::distributed::Triangulation objects and " "associated DoF handler objects, asking for any subdomain other " "than the locally owned one does not make sense."));

	std::vector<types::global_dof_index> dofs_on_this_cell;
	std::vector<types::global_dof_index> dofs_on_other_cell;
	dofs_on_this_cell.reserve(max_dofs_per_cell(dof));
	dofs_on_other_cell.reserve(max_dofs_per_cell(dof));
	typename DH::active_cell_iterator cell = dof.begin_active(), endc =
			dof.end();

	// TODO: in an old implementation, we used user flags before to tag
	// faces that were already touched. this way, we could reduce the work
	// a little bit. now, we instead add only data from one side. this
	// should be OK, but we need to actually verify it.

	// In case we work with a distributed sparsity pattern of Trilinos
	// type, we only have to do the work if the current cell is owned by
	// the calling processor. Otherwise, just continue.
	for (; cell != endc; ++cell)
		if (((subdomain_id == numbers::invalid_subdomain_id)
				|| (subdomain_id == cell->subdomain_id()))
				&& cell->is_locally_owned()) {
			const unsigned int n_dofs_on_this_cell =
					cell->get_fe().dofs_per_cell;
			dofs_on_this_cell.resize(n_dofs_on_this_cell);
			cell->get_dof_indices(dofs_on_this_cell);

			// make sparsity pattern for this cell. if no constraints pattern
			// was given, then the following call acts as if simply no
			// constraints existed
			constraints.add_entries_local_to_global(dofs_on_this_cell, sparsity,
					keep_constrained_dofs);

			for (unsigned int face = 0;
					face < GeometryInfo<DH::dimension>::faces_per_cell;
					++face) {
				typename DH::face_iterator cell_face = cell->face(face);
				if (!cell->at_boundary(face)) {
					typename DH::level_cell_iterator neighbor = cell->neighbor(
							face);

					// in 1d, we do not need to worry whether the neighbor
					// might have children and then loop over those children.
					// rather, we may as well go straight to to cell behind
					// this particular cell's most terminal child
					if (DH::dimension == 1)
						while (neighbor->has_children())
							neighbor = neighbor->child(face == 0 ? 1 : 0);

					if (neighbor->has_children()) {
						for (unsigned int sub_nr = 0;
								sub_nr != cell_face->number_of_children();
								++sub_nr) {
							const typename DH::level_cell_iterator sub_neighbor =
									cell->neighbor_child_on_subface(face,
											sub_nr);

							const unsigned int n_dofs_on_neighbor =
									sub_neighbor->get_fe().dofs_per_cell;
							dofs_on_other_cell.resize(n_dofs_on_neighbor);
							sub_neighbor->get_dof_indices(dofs_on_other_cell);

							constraints.add_entries_local_to_global(
									dofs_on_this_cell, dofs_on_other_cell,
									sparsity, keep_constrained_dofs);
							constraints.add_entries_local_to_global(
									dofs_on_other_cell, dofs_on_this_cell,
									sparsity, keep_constrained_dofs);
							// only need to add this when the neighbor is not
							// owned by the current processor, otherwise we add
							// the entries for the neighbor there
							if (sub_neighbor->subdomain_id()
									!= cell->subdomain_id())
								constraints.add_entries_local_to_global(
										dofs_on_other_cell, sparsity,
										keep_constrained_dofs);
						}
					} else {
						// Refinement edges are taken care of by coarser
						// cells
						if (cell->neighbor_is_coarser(face)
								&& neighbor->subdomain_id()
										== cell->subdomain_id())
							continue;

						const unsigned int n_dofs_on_neighbor =
								neighbor->get_fe().dofs_per_cell;
						dofs_on_other_cell.resize(n_dofs_on_neighbor);

						neighbor->get_dof_indices(dofs_on_other_cell);

						constraints.add_entries_local_to_global(
								dofs_on_this_cell, dofs_on_other_cell, sparsity,
								keep_constrained_dofs);

						// only need to add these in case the neighbor cell
						// is not locally owned - otherwise, we touch each
						// face twice and hence put the indices the other way
						// around
						if (!cell->neighbor(face)->active()
								|| (neighbor->subdomain_id()
										!= cell->subdomain_id())) {
							constraints.add_entries_local_to_global(
									dofs_on_other_cell, dofs_on_this_cell,
									sparsity, keep_constrained_dofs);
							if (neighbor->subdomain_id()
									!= cell->subdomain_id())
								constraints.add_entries_local_to_global(
										dofs_on_other_cell, sparsity,
										keep_constrained_dofs);
						}
					}
				}
			}
		}
}

template<class DH, class SparsityPattern>
void make_sparser_flux_sparsity_pattern(const DH &dof,
		SparsityPattern &sparsity) {
	ConstraintMatrix constraints;
	make_sparser_flux_sparsity_pattern(dof, sparsity, constraints);
}

} /* namepace DealIIExtensions */

} /* namespace natrium */

typedef CompressedSparsityPattern SP;
//for (SP : SPARSITY_PATTERNS; deal_II_dimension : DIMENSIONS)
//for (size_t deal_II_dimension = 1; deal_II_dimension < 4; deal_II_dimension++)
//{
template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<2>, SP>(
		const DoFHandler<2> &dof, SP &sparsity);

template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<2>, SP>(
		const DoFHandler<2> &dof, SP &sparsity,
		const ConstraintMatrix &constraints, const bool, const unsigned int);

template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<3>, SP>(
		const DoFHandler<3> &dof, SP &sparsity);

template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<3>, SP>(
		const DoFHandler<3> &dof, SP &sparsity,
		const ConstraintMatrix &constraints, const bool, const unsigned int);

//}

