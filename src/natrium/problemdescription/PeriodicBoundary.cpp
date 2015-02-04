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

template<size_t dim> PeriodicBoundary<dim>::PeriodicBoundary(
		size_t boundaryIndicator1, size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<dim> > triangulation) :
		m_boundaryIndicator1(boundaryIndicator1), m_boundaryIndicator2(
				boundaryIndicator2) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2) {
		throw PeriodicBoundaryNotPossible(
				"The boundary indicators defining the periodic boundaries must not be equal to each other.");
	}

	m_triangulation = triangulation;

	if (dim == 2) {
		// create vertex points
		dealii::Point<2> beginLine1(0.0, 0.0);
		dealii::Point<2> endLine1(0.0, 0.0);
		dealii::Point<2> beginLine2(0.0, 0.0);
		dealii::Point<2> endLine2(0.0, 0.0);

		// calculate the positions of the vertex points(
		std::string errorMessage1;
		bool areLines = BoundaryTools::getInterfacialLinesByBoundaryIndicator(
				boundaryIndicator1, boundaryIndicator2, triangulation,
				beginLine1, endLine1, beginLine2, endLine2, errorMessage1);
		if (not areLines)
			throw PeriodicBoundaryNotPossible(errorMessage1);

		m_beginLine1 = beginLine1;
		m_beginLine2 = beginLine2;
		m_endLine1 = endLine1;
		m_endLine2 = endLine2;

		// check if positions of the interfaces are OK;
		// else throw PeriodicBoundaryNotPossible
		std::string errorMessage2;
		bool isParallel = BoundaryTools::checkParallelLines(m_beginLine1,
				m_endLine1, m_beginLine2, m_endLine2, errorMessage2);
		if (not isParallel)
			throw PeriodicBoundaryNotPossible(errorMessage2);
	}
} /* Constructor 2 */
// The template Parameter has to be made explicit in order for the code to compile
template PeriodicBoundary<2>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation);

template<size_t dim> bool PeriodicBoundary<dim>::isFaceInBoundary(
		const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
		size_t faceBoundaryIndicator) const {
	// first condition: cell map has a key <cell>
	if (m_cells.count(cell) == 0) {
		return false;
	}
	// second condition: the face has the right boundary indicator
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

template<> void PeriodicBoundary<2>::createCellMap(
		const dealii::DoFHandler<2>& doFHandler) {

	// Make iterators over active faces
	dealii::DoFHandler<2>::active_cell_iterator currentCell =
			doFHandler.begin_active();
	dealii::DoFHandler<2>::active_cell_iterator lastCell = doFHandler.end();

	// The key of these cells is the distance from the begin point of the respective line.
	// The second element of the value pair is the local face id of the face which belongs to the boundary.
	std::map<double,
			std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t>, own_double_less > cellsAtBoundary1;
	std::map<double, std::pair<dealii::DoFHandler<2>::cell_iterator, size_t>, own_double_less > cellsAtBoundary2;

	// iterate over all active cells and sort them
	for (; currentCell != lastCell; ++currentCell) {
		if (currentCell->at_boundary()) {
			for (size_t i = 0; i < dealii::GeometryInfo<2>::faces_per_cell;
					i++) {
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator1) {
					double key = (currentCell->face(i))->center().distance(
							m_beginLine1);
					cellsAtBoundary1.insert(
							std::make_pair(key,
									std::make_pair(currentCell,
											dealii::GeometryInfo<2>::opposite_face[i])));
				}
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator2) {
					double key = (currentCell->face(i))->center().distance(
							m_beginLine2);
					cellsAtBoundary2.insert(
							std::make_pair(key,
									std::make_pair(currentCell,
											dealii::GeometryInfo<2>::opposite_face[i])));
				}
			}
		}
	}

	// make sure that both boundaries have the same number of adjacent cells
	if (cellsAtBoundary1.size() != cellsAtBoundary2.size()) {
		throw PeriodicBoundaryNotPossible(
				"The boundaries do not have the same number of adjacent cells.");
	}

	// store the cell ids and their respective opposite cell in map
	m_cells.clear();
	std::map<double,
			std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t> >::iterator atBoundary1 =
			cellsAtBoundary1.begin();
	std::map<double, std::pair<dealii::DoFHandler<2>::cell_iterator, size_t> >::iterator atBoundary2 =
			cellsAtBoundary2.begin();
	for (; atBoundary1 != cellsAtBoundary1.end(); atBoundary1++) {
		// assert that the face discretizations on both lines are equal
		if (cellsAtBoundary2.count(atBoundary1->first) == 0) {
			throw PeriodicBoundaryNotPossible(
					"The discretizations of opposite periodic boundaries do not coincide. This version of the NATriuM solver does only work with equal discretizations.");
		}
		// add to cells
		m_cells.insert(
				std::make_pair(atBoundary1->second.first, atBoundary2->second));
		m_cells.insert(
				std::make_pair(atBoundary2->second.first, atBoundary1->second));
		atBoundary2++;
	}

} /* createMap */
template<> void PeriodicBoundary<3>::createCellMap(
		const dealii::DoFHandler<3>& doFHandler) {
}

template<size_t dim> size_t PeriodicBoundary<dim>::getOppositeCellAtPeriodicBoundary(
		const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
		typename dealii::DoFHandler<dim>::cell_iterator & neighborCell) const {

	if (m_cells.size() == 0) {
		throw PeriodicBoundaryNotPossible(
				"CreateMap has to be called before getOppositeCellAtPeriodicBoundary.");
	}
// assert that the given cell is at the boundary
	if (m_cells.count(cell) == 0) {
		throw PeriodicBoundaryNotPossible(
				"The cell does not belong to the boundary.");
	}

	neighborCell = m_cells.at(cell).first;
	return m_cells.at(neighborCell).second;
}
// The template parameter has to be made explicit in order for the code to compile
template size_t PeriodicBoundary<2>::getOppositeCellAtPeriodicBoundary(
		const dealii::DoFHandler<2>::active_cell_iterator & cell,
		dealii::DoFHandler<2>::cell_iterator & neighborCell) const;
template size_t PeriodicBoundary<3>::getOppositeCellAtPeriodicBoundary(
		const dealii::DoFHandler<3>::active_cell_iterator & cell,
		dealii::DoFHandler<3>::cell_iterator & neighborCell) const;

template<size_t dim> void PeriodicBoundary<dim>::addToSparsityPattern(
		dealii::BlockCompressedSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const {

	// THIS FUNCTION IS NOT USED!!! See DealIIExtensions module for details



	// ConstraintMatrix can be used for a more efficient distribution to global sparsity patterns
	const dealii::ConstraintMatrix constraints;

// add periodic boundaries to intermediate flux sparsity pattern
	vector<dealii::types::global_dof_index> doFIndicesAtCell1(dofs_per_cell);
	vector<dealii::types::global_dof_index> doFIndicesAtCell2(dofs_per_cell);
	// for all blocks (it is important to have this loop in the outer part
	// because otherwise it is very Cache-inefficient)
	for (size_t I = 0; I < n_blocks; I++) {
		//minimize calls of block(), because expensive
		dealii::CompressedSparsityPattern& block = cSparse.block(I, I);
		typename std::map<
				typename dealii::DoFHandler<dim>::active_cell_iterator,
				std::pair<typename dealii::DoFHandler<dim>::cell_iterator,
						size_t> >::const_iterator element = m_cells.begin();
		// for each cells belonging to the periodic boundary
		for (; element != m_cells.end(); element++) {
			element->first->get_dof_indices(doFIndicesAtCell1);
			element->second.first->get_dof_indices(doFIndicesAtCell2);
			// couple all dofs at boundary 1 with dofs at boundary 2
			constraints.add_entries_local_to_global(doFIndicesAtCell1, doFIndicesAtCell2, block, true);
			// TODO only couple the ones which are nonzero at the face
			// TODO remove the INVARIANT "discretization at boundary 1 = discretization at boundary 2"
			//      e.g. by mapping, allowing more than one periodic neighbor, ...
			// iterate over rows
			//for (size_t j = 0; j < dofs_per_cell; j++) {
			//	block.add_entries (doFIndicesAtCell1[j], doFIndicesAtCell2.begin(), doFIndicesAtCell2.end());
			//}
		}
	}


}
template void PeriodicBoundary<2>::addToSparsityPattern(
		dealii::BlockCompressedSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const;
template void PeriodicBoundary<3>::addToSparsityPattern(
		dealii::BlockCompressedSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const;

} /* namespace natrium */
