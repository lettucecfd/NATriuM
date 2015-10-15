/*
 * DealiiExtensions.h
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 */

#ifndef DEALIIEXTENSIONS_H_
#define DEALIIEXTENSIONS_H_

#include <vector>
#include <set>
#include <map>

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_iterator.h>

#include "../problemdescription/PeriodicBoundary.h"

namespace natrium {

template<size_t dim, size_t spacedim = dim>
using FacePair = dealii::GridTools::PeriodicFacePair< dealii::TriaIterator<dealii::CellAccessor<dim, spacedim> > >;


template<size_t dim, size_t spacedim = dim>
using PeriodicCellMap = std::map<
			dealii::TriaIterator< dealii::CellAccessor<dim, spacedim> >,
			FacePair<dim, spacedim>
	>;

// forward declarations
class PeriodicBoundaryNotPossible;
template <size_t dim>
class BoundaryCollection;

namespace DealIIExtensions {

//DEAL_II_NAMESPACE_OPEN

using namespace dealii;


/**
 * @short Like dealii::DoFTools::make_flux_sparsity_pattern but does only create
 * 		  non-zero entries for the DoFs situated on faces
 */
template<class DH, class SparsityPattern>
void
make_sparser_flux_sparsity_pattern(const DH &dof, SparsityPattern &sparsity,
		const ConstraintMatrix &constraints,
		const BoundaryCollection<DH::dimension>& boundaries =
				natrium::BoundaryCollection<DH::dimension>(),
		FEFaceValues<DH::dimension>* fe_face = NULL,
		const bool keep_constrained_dofs = true,
		const types::subdomain_id subdomain_id = numbers::invalid_unsigned_int);

template<class DH, class SparsityPattern>
void
make_sparser_flux_sparsity_pattern(const DH &dof, SparsityPattern &sparsity,
		const natrium::BoundaryCollection<DH::dimension>& boundaries =
				BoundaryCollection<DH::dimension>(),
		FEFaceValues<DH::dimension>* fe_face = NULL);

template<typename DH>
void make_periodicity_map_dg(const typename DH::cell_iterator &cell_1,
		const typename identity<typename DH::cell_iterator>::type &cell_2, size_t face_nr_1,
		size_t face_nr_2, PeriodicCellMap<DH::dimension>& cell_map,
		const bool face_orientation, const bool face_flip,
		const bool face_rotation);

template<typename DH>
void make_periodicity_map_dg(
		const std::vector<
				GridTools::PeriodicFacePair<typename DH::cell_iterator> > &periodic_faces,
		PeriodicCellMap<DH::dimension>& cell_map);

template<typename DH>
void make_periodicity_map_dg(const DH &dof_handler,
		size_t b_id1, size_t b_id2,
		const int direction, PeriodicCellMap<DH::dimension>& cell_map);

//DEAL_II_NAMESPACE_CLOSE

} /* namespace DealIIExtensions */

} /* namespace natrium */

#endif /* DEALIIEXTENSIONS_H_ */

