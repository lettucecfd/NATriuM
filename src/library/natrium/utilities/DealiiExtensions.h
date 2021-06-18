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

#include <deal.II/grid/tria.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_iterator.h>

//#include "../boundaries/PeriodicBoundary.h"

namespace natrium {

template<size_t dim, size_t spacedim = dim>
using FacePair = dealii::GridTools::PeriodicFacePair< dealii::TriaIterator<dealii::CellAccessor<dim, spacedim> > >;

/**
 * @short	A std::map that maps cells to face pairs
 * @note	The order of cells in the Face Pair is of no meaning
 */
template<size_t dim, size_t spacedim = dim>
using PeriodicCellMap = std::map<
dealii::TriaIterator< dealii::CellAccessor<dim, spacedim> >,
FacePair<dim, spacedim>
>;

// forward declarations
class PeriodicBoundaryNotPossible;
template<size_t dim>
class BoundaryCollection;

/** @short
 * some extensions to deal.ii
 */
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
		const AffineConstraints<double> &constraints,
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

/**
 *  @short 	Gathers cell pairs at a periodic boundary. This function starts at the coarsest level and
 *  		recursively visits the subcells of boundary cells. At the active level, cell pairs are added
 *  		to the cell_map.
 *  @param[in] cell_1 		cell at the one side of a periodic boundary
 *  @param[in] cell_2 		cell at the other side of a periodic boundary
 *  @param[in] face_nr_1	boundary face number of the cell_1
 *  @param[in] face_nr_2	boundary face number of the cell_2
 *  @param[out] cell_map	A map that stores cells and faces at the periodic boundary
 *  @param[in] face_orientation		see deal.ii Glossary
 *  @param[in] face_flip			see deal.ii Glossary
 *  @param[in] face_rotation		see deal.ii Glossary
 *  @note The implementation of this function is similar to dealii::make_periodicity_constraints
 */
template<typename DH>
void make_periodicity_map_dg(const typename DH::cell_iterator &cell_1,
		const typename identity<typename DH::cell_iterator>::type &cell_2,
		size_t face_nr_1, size_t face_nr_2,
		PeriodicCellMap<DH::dimension>& cell_map, const bool face_orientation,
		const bool face_flip, const bool face_rotation);

/**
 * @short	High-level version of the first function, starting from PeriodicFacePairs
 * @param[in]	periodic_faces 	a vector of PeriodicFacePairs. They can be obtained by
 * 								calling collect_periodic_faces.
 * @param[out]	cell_map		A map that stores cells and faces at the periodic boundary
 */
template<typename DH>
void make_periodicity_map_dg(
		const std::vector<
				GridTools::PeriodicFacePair<typename DH::cell_iterator> > &periodic_faces,
		PeriodicCellMap<DH::dimension>& cell_map);

/**
 * @short	Another high-level version of the first function, starting from boundary indicators
 * @param[in]	dof_handler 	a DoFHandler object
 * @param[in]	b_id1			boundary id at the first side of the periodic boundary
 * @param[in]	b_id2			boundary id at the second side of the periodic boundary
 * @param[out]	cell_map		A map that stores cells and faces at the periodic boundary
 * @note 	Creating periodic boundaries requires great care in the order of operations (at least for
 * 			Meshes of type parallel::distributed::Triangulation<dim>)
 * 			-# create unrefined mesh
 * 			-# call collect_periodic_faces<Mesh> to find out the periodic couplings at the coarsest level
 * 			-# call add_periodicity to create a layer of ghost nodes across the periodic boundary
 * 			-# call this function (make_periodicity_map_dg) to recursively find out the periodic couplings
 * 			   on the active levels
 */
template<typename DH>
void make_periodicity_map_dg(const DH &dof_handler, size_t b_id1, size_t b_id2,
		const int direction, PeriodicCellMap<DH::dimension>& cell_map);

/**
 * @short This function extracts those degrees of freedom whose shape functions
 * are nonzero on at least part of the selected boundary.
 * For continuous elements, this is exactly the set of shape functions whose
 * degrees of freedom are defined on boundary faces. On the other hand, if the
 * finite element in used is a discontinuous element, all degrees of freedom
 * are defined in the inside of cells and consequently none would be boundary
 * degrees of freedom. Several of those would have shape functions that are
 * nonzero on the boundary, however. This function therefore extracts all those
 * for which the FiniteElement::has_support_on_face function says that it is
 * nonzero on any face on one of the selected boundary parts.
 *
 * @note This function is the same as the deal.ii function with the same name;
 * except for one line which restricts the function to the locally owned cells.
 */
template<class DH>
void extract_dofs_with_support_on_boundary(const DH &dof_handler,
		const ComponentMask &component_mask, std::vector<bool> &selected_dofs,
		const std::set<types::boundary_id> &boundary_ids);

/**
 * @short Set boundary ids for all faces that have their vertices on the hyperplane x[component]=value
 * @param mesh the mesh
 * @param component 0 for x, 1 for y, 2 for z
 * @param value the component value
 * @param boundary_id the boundary id that is assigned to the faces that are in the hyperpalen
 * @param tol tolerance for comparison with value
 *
 */
template<size_t dim>
void set_boundary_ids_at_hyperplane(dealii::Triangulation<dim>& mesh,
		 unsigned int component, double value, size_t boundary_id, double tol = 1e-10);


//DEAL_II_NAMESPACE_CLOSE

} /* namespace DealIIExtensions */

} /* namespace natrium */



#endif /* DEALIIEXTENSIONS_H_ */

