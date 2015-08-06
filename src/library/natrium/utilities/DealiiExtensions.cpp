/*
 * DealiiExtensions.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 */

#include "DealiiExtensions.h"

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

#include <algorithm>
#include <numeric>

using namespace dealii;
using namespace DoFTools;

namespace natrium {

namespace DealIIExtensions {

template<class DH, class SparsityPattern>
void make_sparser_flux_sparsity_pattern(const DH &dof,
		SparsityPattern &sparsity, const ConstraintMatrix &constraints,
		const natrium::BoundaryCollection<DH::dimension>& boundaries,
		FEFaceValues<DH::dimension>* fe_face, const bool keep_constrained_dofs,
		const types::subdomain_id subdomain_id)

		{

	const types::global_dof_index n_dofs = dof.n_dofs();

	AssertDimension(sparsity.n_rows(), n_dofs);
	AssertDimension(sparsity.n_cols(), n_dofs);
	if (fe_face != NULL) {
		Assert(fe_face->get_fe().has_support_points(),
				ExcMessage ("Sparser flux sparsity pattern makes only sense for elements with support points"));
	}

	// If we have a distributed::Mesh only allow locally_owned
	// subdomain. Not setting a subdomain is also okay, because we skip
	// ghost cells in the loop below.
	Assert(
			(dof.get_tria().locally_owned_subdomain() == numbers::invalid_subdomain_id) || (subdomain_id == numbers::invalid_subdomain_id) || (subdomain_id == dof.get_tria().locally_owned_subdomain()),
			ExcMessage ("For parallel::distributed::Mesh objects and " "associated DoF handler objects, asking for any subdomain other " "than the locally owned one does not make sense."));

	std::vector<types::global_dof_index> dofs_on_this_cell;
	std::vector<types::global_dof_index> dofs_on_other_cell;
	std::vector<types::global_dof_index> dofs_on_this_face;
	std::vector<types::global_dof_index> dofs_on_other_face;
	dofs_on_this_cell.reserve(max_dofs_per_cell(dof));
	dofs_on_other_cell.reserve(max_dofs_per_cell(dof));
	dofs_on_this_face.reserve(max_dofs_per_cell(dof));
	dofs_on_other_face.reserve(max_dofs_per_cell(dof));
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
				typename DH::face_iterator this_face = cell->face(face);
				typename DH::face_iterator other_face;
				typename DH::cell_iterator neighbor;
				typename DH::cell_iterator test_cell;
				unsigned int neighbor_face = 1000;
				unsigned int test_face;
				if (cell->at_boundary(face)) {
					// check if cell is in periodic boundary
					for (typename BoundaryCollection<DH::dimension>::ConstPeriodicIterator periodic =
							boundaries.getPeriodicBoundaries().begin();
							periodic != boundaries.getPeriodicBoundaries().end();
							periodic++) {
						if (periodic->second->isFaceInBoundary(cell,
								this_face->boundary_indicator())) {
							neighbor_face =
									periodic->second->getOppositeCellAtPeriodicBoundary(
											cell, neighbor);
							// if the face is not really at the boundary: do not consider cell
							/*test_face =
									periodic->second->getOppositeCellAtPeriodicBoundary(
											neighbor, test_cell);
							if (not ((test_cell == cell) && (test_face == face))) {
								neighbor_face = 1000;
								break;
							} else {*/
								other_face = neighbor->face(neighbor_face);
								break;
							//}
						}
					}
					if (neighbor_face == 1000) {
						continue;
					}
				} else {
					neighbor = cell->neighbor(face);
					neighbor_face = 1000; // indicator that it has to be assembled, later
				}

				// specify whether pairwise coupling is valid
				bool pairwise_coupling_valid = (fe_face != NULL);

				dofs_on_this_face.clear();
				dofs_on_other_face.clear();
				if (!pairwise_coupling_valid) {
					// get face dofs
					for (size_t i = 0; i < n_dofs_on_this_cell; i++) {
						if (cell->get_fe().has_support_on_face(i, face)) {
							dofs_on_this_face.push_back(
									dofs_on_this_cell.at(i));
						}
					}
				}
				// in 1d, we do not need to worry whether the neighbor
				// might have children and then loop over those children.
				// rather, we may as well go straight to to cell behind
				// this particular cell's most terminal child
				if (DH::dimension == 1)
					while (neighbor->has_children())
						neighbor = neighbor->child(face == 0 ? 1 : 0);

				if (neighbor->has_children()) {
					for (unsigned int sub_nr = 0;
							sub_nr != this_face->number_of_children();
							++sub_nr) {
						const typename DH::cell_iterator sub_neighbor =
								cell->neighbor_child_on_subface(face, sub_nr);

						const unsigned int n_dofs_on_neighbor =
								sub_neighbor->get_fe().dofs_per_cell;
						dofs_on_other_cell.resize(n_dofs_on_neighbor);
						sub_neighbor->get_dof_indices(dofs_on_other_cell);

						// identify which sub_neighbor face is child to this face
						unsigned int sub_neighbor_face;
						for (sub_neighbor_face = 0;
								sub_neighbor_face
										< GeometryInfo<DH::dimension>::faces_per_cell;
								++sub_neighbor_face) {
							other_face = sub_neighbor->face(sub_neighbor_face);
							if (sub_neighbor->neighbor(sub_neighbor_face)
									== cell) {
								break;
							}
							Assert(
									sub_neighbor_face + 1 < GeometryInfo<DH::dimension>::faces_per_cell,
									ExcMessage ("Neighbor face was not found, but needed for constructing the sparsity pattern"));
						}

						// Couple all dofs on common face
						dofs_on_other_face.clear();
						for (size_t i = 0; i < n_dofs_on_neighbor; i++) {
							if (sub_neighbor->get_fe().has_support_on_face(i,
									sub_neighbor_face)) {
								dofs_on_other_face.push_back(
										dofs_on_other_cell.at(i));
							}
						}
						Assert(
								dofs_on_this_face.size() * dofs_on_other_face.size() > 0,
								ExcMessage ("Size of at least one dof vector is 0."));

						// Add entries to sparsity pattern
						constraints.add_entries_local_to_global(
								dofs_on_this_face, dofs_on_other_face, sparsity,
								keep_constrained_dofs);
						constraints.add_entries_local_to_global(
								dofs_on_other_face, dofs_on_this_face, sparsity,
								keep_constrained_dofs);
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
					if (!cell->at_boundary(face)) {
						if (cell->neighbor_is_coarser(face)
								&& neighbor->subdomain_id()
										== cell->subdomain_id())
							continue;
					}
					const unsigned int n_dofs_on_neighbor =
							neighbor->get_fe().dofs_per_cell;
					dofs_on_other_cell.resize(n_dofs_on_neighbor);
					neighbor->get_dof_indices(dofs_on_other_cell);

					// identify which neighbor face belongs to this face
					if (neighbor_face == 1000) {
						for (neighbor_face = 0;
								neighbor_face
										< GeometryInfo<DH::dimension>::faces_per_cell;
								++neighbor_face) {
							other_face = neighbor->face(neighbor_face);
							if (*other_face == *this_face) {
								break;
							}
							Assert(
									neighbor_face + 1 < GeometryInfo<DH::dimension>::faces_per_cell,
									ExcMessage ("Neighbor face was not found, but needed for constructing the sparsity pattern"));
						}
					}

					if (!pairwise_coupling_valid) {
						// Method 1) Couple all dofs on common face
						dofs_on_other_face.clear();
						for (size_t i = 0; i < n_dofs_on_neighbor; i++) {
							if (neighbor->get_fe().has_support_on_face(i,
									neighbor_face)) {
								dofs_on_other_face.push_back(
										dofs_on_other_cell.at(i));
							}
						}
						Assert(
								dofs_on_this_face.size() * dofs_on_other_face.size() > 0,
								ExcMessage ("Size of at least one dof vector is 0."));

						// Add entries to sparsity pattern
						constraints.add_entries_local_to_global(
								dofs_on_this_face, dofs_on_other_face, sparsity,
								keep_constrained_dofs);
					} else {
						// Method 2) if possible: make unique relation between neighboring dofs
						fe_face->reinit(cell, face);
						// bring dofs at this face into right order
						for (size_t i = 0; i < n_dofs_on_this_cell; i++) {
							int unique = 0;
							for (size_t q = 0; q < fe_face->n_quadrature_points;
									q++) {
								if (fe_face->shape_value(i, q) > 1e-10) {
									unique += 1;
									dofs_on_this_face.push_back(
											dofs_on_this_cell.at(i));
								}
							}
							// Test, if the relationship doF <-> quadrature points in unique
							assert(unique <= 1);
						}
						// bring dofs at other face in right order
						fe_face->reinit(neighbor, neighbor_face);
						for (size_t i = 0; i < n_dofs_on_neighbor; i++) {
							int unique = 0;
							for (size_t q = 0; q < fe_face->n_quadrature_points;
									q++) {
								if (fe_face->shape_value(i, q) > 1e-10) {
									unique += 1;
									dofs_on_other_face.push_back(
											dofs_on_other_cell.at(i));
								}
							}
							// Test, if the relationship doF <-> quadrature points in unique
							assert(unique <= 1);
						}

						if (2 == DH::dimension){
						AssertDimension(dofs_on_this_face.size(),
								sqrt(n_dofs_on_this_cell));
						} else if (3 == DH::dimension){
							AssertDimension(dofs_on_this_face.size()* sqrt(dofs_on_this_face.size()),
									n_dofs_on_this_cell);
						}
						AssertDimension(dofs_on_other_face.size(),
								dofs_on_this_face.size());

						// couple only individual dofs with each other
						std::vector<types::global_dof_index> dof_this(1);
						std::vector<types::global_dof_index> dof_other(1);
						for (size_t i = 0; i < dofs_on_this_face.size(); i++) {
							dof_this.at(0) = dofs_on_this_face.at(i);
							dof_other.at(0) = dofs_on_other_face.at(i);

							// Add entries to sparsity pattern
							constraints.add_entries_local_to_global(dof_this,
									dof_other, sparsity, keep_constrained_dofs);
						}
					}
					// only need to add these in case the neighbor cell
					// is not locally owned - otherwise, we touch each
					// face twice and hence put the indices the other way
					// around
					if (!(neighbor->active())
							|| (neighbor->subdomain_id() != cell->subdomain_id())) {
						if (!pairwise_coupling_valid) {
							constraints.add_entries_local_to_global(
									dofs_on_other_face, dofs_on_this_face,
									sparsity, keep_constrained_dofs);
						} else {
							// couple only individual dofs with each other
							std::vector<types::global_dof_index> dof_this(1);
							std::vector<types::global_dof_index> dof_other(1);
							for (size_t i = 0; i < dofs_on_this_face.size();
									i++) {
								dof_this.at(0) = dofs_on_this_face.at(i);
								dof_other.at(0) = dofs_on_other_face.at(i);

								// Add entries to sparsity pattern
								constraints.add_entries_local_to_global(
										dof_this, dof_other, sparsity,
										keep_constrained_dofs);
							}
						}
						if (neighbor->subdomain_id() != cell->subdomain_id())
							constraints.add_entries_local_to_global(
									dofs_on_other_cell, sparsity,
									keep_constrained_dofs);
					}
				}
			}
		}
}

template<class DH, class SparsityPattern>
void make_sparser_flux_sparsity_pattern(const DH &dof,
		SparsityPattern &sparsity,
		const natrium::BoundaryCollection<DH::dimension>& boundaries,
		FEFaceValues<DH::dimension>* fe_face) {
	ConstraintMatrix constraints;
	make_sparser_flux_sparsity_pattern(dof, sparsity, constraints, boundaries,
			fe_face);
}

} /* namepace DealIIExtensions */

} /* namespace natrium */

typedef DynamicSparsityPattern SP;
//for (SP : SPARSITY_PATTERNS; deal_II_dimension : DIMENSIONS)
//for (size_t deal_II_dimension = 1; deal_II_dimension < 4; deal_II_dimension++)
//{
template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<2>, SP>(
		const DoFHandler<2> &dof, SP &sparsity,
		const natrium::BoundaryCollection<2>& boundaries,
		FEFaceValues<2>* fe_face);

template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<2>, SP>(
		const DoFHandler<2> &dof, SP &sparsity,
		const ConstraintMatrix &constraints,
		const natrium::BoundaryCollection<2>& boundaries,
		FEFaceValues<2>* fe_face, const bool, const unsigned int);

template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<3>, SP>(
		const DoFHandler<3> &dof, SP &sparsity,
		const natrium::BoundaryCollection<3>& boundaries,
		FEFaceValues<3>* fe_face);

template void
natrium::DealIIExtensions::make_sparser_flux_sparsity_pattern<DoFHandler<3>, SP>(
		const DoFHandler<3> &dof, SP &sparsity,
		const ConstraintMatrix &constraints,
		const natrium::BoundaryCollection<3>& boundaries,
		FEFaceValues<3>* fe_face, const bool, const unsigned int);

//}

