/**
 * @file PeriodicBoundary.cpp
 * @short Description of a periodic boundary
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PeriodicBoundary.h"

#include <iterator>

#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/lac/constraint_matrix.h"

#include "../utilities/BasicNames.h"
#include "../utilities/Math.h"

#include "BoundaryTools.h"

namespace natrium {

// The template Parameter has to be made explicit in order for the code to compile
template<size_t dim> PeriodicBoundary<dim>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2, size_t direction,
		shared_ptr<Mesh<dim> > triangulation) :
		m_boundaryIndicator1(boundaryIndicator1), m_boundaryIndicator2(
				boundaryIndicator2), m_direction(direction) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2) {
		throw PeriodicBoundaryNotPossible(
				"The boundary indicators defining the periodic boundaries must not be equal to each other.");
	}

	m_triangulation = triangulation;

} /* Constructor 2 */
template PeriodicBoundary<2>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2, size_t direction,
		shared_ptr<Mesh<2> > triangulation);
template PeriodicBoundary<3>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2, size_t direction,
		shared_ptr<Mesh<3> > triangulation);

template<size_t dim> bool PeriodicBoundary<dim>::isFaceInBoundary(
		const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
		size_t faceBoundaryIndicator) const {
	/*
	// first condition: cell map has a key <cell>
	if (m_cells.count(cell) == 0) {
		return false;
	}
	// second condition: the face has the right boundary indicator
	 */
	assert(m_cells.count(cell) != 0);
	if (faceBoundaryIndicator == m_boundaryIndicator1) {
		return true;
	}
	if (faceBoundaryIndicator == m_boundaryIndicator2) {
		return true;
	}
	return false;
}
// The template Parameter has to be made explicit in order for the code to compile
template bool PeriodicBoundary<2>::isFaceInBoundary(
		const typename dealii::DoFHandler<2>::active_cell_iterator & cell,
		size_t faceBoundaryIndicator) const;
template bool PeriodicBoundary<3>::isFaceInBoundary(
		const typename dealii::DoFHandler<3>::active_cell_iterator & cell,
		size_t faceBoundaryIndicator) const;

template<size_t dim> PeriodicBoundary<dim>::~PeriodicBoundary() {
}
// The template Parameter has to be made explicit in order for the code to compile
template PeriodicBoundary<2>::~PeriodicBoundary();
template PeriodicBoundary<3>::~PeriodicBoundary();

template<size_t dim> void PeriodicBoundary<dim>::createCellMap(
		const dealii::DoFHandler<dim>& doFHandler) {

	DealIIExtensions::make_periodicity_map_dg(doFHandler,
			m_boundaryIndicator1, m_boundaryIndicator2, m_direction, m_cells);

} /* createMap */
template void PeriodicBoundary<2>::createCellMap(
		const dealii::DoFHandler<2>& doFHandler);
template void PeriodicBoundary<3>::createCellMap(
		const dealii::DoFHandler<3>& doFHandler);


template<size_t dim> size_t PeriodicBoundary<dim>::getOppositeCellAtPeriodicBoundary(
		const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
		typename dealii::DoFHandler<dim>::cell_iterator & neighborCell) const {

	if (m_cells.size() == 0) {
		throw PeriodicBoundaryNotPossible(
				"CreateCellMap has to be called before getOppositeCellAtPeriodicBoundary.");
	}
// assert that the given cell is at the boundary
	if (m_cells.count(cell) == 0) {
		throw PeriodicBoundaryNotPossible(
				"The cell does not belong to the boundary.");
	}
	FacePair<dim> face_pair = m_cells.at(cell);
	if (face_pair.cell[0] == cell){
		neighborCell = face_pair.cell[1];
		return face_pair.face_idx[1];
	} else {
		assert (face_pair.cell[1] == cell);
		neighborCell = face_pair.cell[0];
		return face_pair.face_idx[0];
	}
}
// The template parameter has to be made explicit in order for the code to compile
template size_t PeriodicBoundary<2>::getOppositeCellAtPeriodicBoundary(
		const dealii::DoFHandler<2>::active_cell_iterator & cell,
		dealii::DoFHandler<2>::cell_iterator & neighborCell) const;
template size_t PeriodicBoundary<3>::getOppositeCellAtPeriodicBoundary(
		const dealii::DoFHandler<3>::active_cell_iterator & cell,
		dealii::DoFHandler<3>::cell_iterator & neighborCell) const;

template<size_t dim> void PeriodicBoundary<dim>::addToSparsityPattern(
		dealii::BlockDynamicSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const {

	// THIS FUNCTION IS NOT USED!!! See DealIIExtensions module for details
	assert (false);
}
template void PeriodicBoundary<2>::addToSparsityPattern(
		dealii::BlockDynamicSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const;
template void PeriodicBoundary<3>::addToSparsityPattern(
		dealii::BlockDynamicSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const;

} /* namespace natrium */
