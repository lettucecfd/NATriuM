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
template<> PeriodicBoundary<2>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<Mesh<2> > triangulation) :
		m_boundaryIndicator1(boundaryIndicator1), m_boundaryIndicator2(
				boundaryIndicator2) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2) {
		throw PeriodicBoundaryNotPossible(
				"The boundary indicators defining the periodic boundaries must not be equal to each other.");
	}

	m_triangulation = triangulation;

	// create vertex points
	dealii::Point<2> beginLine1(0.0, 0.0);
	dealii::Point<2> endLine1(0.0, 0.0);
	dealii::Point<2> beginLine2(0.0, 0.0);
	dealii::Point<2> endLine2(0.0, 0.0);

	// calculate the positions of the vertex points(
	std::string errorMessage1;
	bool areLines = BoundaryTools::getInterfacialLinesByBoundaryIndicator(
			boundaryIndicator1, boundaryIndicator2, triangulation, beginLine1,
			endLine1, beginLine2, endLine2, errorMessage1);
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

} /* Constructor 2 */
template<> PeriodicBoundary<3>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<Mesh<3> > triangulation) :
		m_boundaryIndicator1(boundaryIndicator1), m_boundaryIndicator2(
				boundaryIndicator2) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2) {
		throw PeriodicBoundaryNotPossible(
				"The boundary indicators defining the periodic boundaries must not be equal to each other.");
	}

	// TODO Meaningfully check if the boundaries are parallel

	m_triangulation = triangulation;
}

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
			std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t>,
			own_double_less> cellsAtBoundary1;
	std::map<double, std::pair<dealii::DoFHandler<2>::cell_iterator, size_t>,
			own_double_less> cellsAtBoundary2;

	// iterate over all active cells and sort them
	for (; currentCell != lastCell; ++currentCell) {
		if (currentCell->at_boundary()) {
			for (size_t i = 0; i < dealii::GeometryInfo<2>::faces_per_cell;
					i++) {
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator1) {
					double key = (currentCell->face(i))->center().distance(
							m_beginLine1);
					if (cellsAtBoundary1.find(key) != cellsAtBoundary1.end()) {
						std::stringstream s;
						s
								<< "Minimum vertex distance < 1e-6; but the periodic boundary detection "
										"tests for cells distances < 1e-6 in own_double_less. That might cause "
										"serious problems." << endl;
						throw PeriodicBoundaryNotPossible(
								"Error in Periodic Boundary: Cells on opposite boundaries not unique."
										"The cell diameter might be < 1e-6, which is not allowed by the function own_double_less.",
								s);
					}
					cellsAtBoundary1.insert(
							std::make_pair(key,
									std::make_pair(currentCell,
											dealii::GeometryInfo<2>::opposite_face[i])));
				}
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator2) {
					double key = (currentCell->face(i))->center().distance(
							m_beginLine2);

					if (cellsAtBoundary2.find(key) != cellsAtBoundary2.end()) {
						std::stringstream s;
						s
								<< "Vertex distance < 1e-6; but the periodic boundary detection "
										"tests for cells distances < 1e-6 in own_double_less. That might cause "
										"serious problems." << endl;
						throw PeriodicBoundaryNotPossible(
								"Error in Periodic Boundary: Cells on opposite boundaries not unique."
										"The cell diameter might be < 1e-6, which is not allowed by the function own_double_less.",
								s);
					}
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
			// Generate lengthy error message
			std::stringstream info;
			info << "Face centers for boundaries " << m_boundaryIndicator1
					<< " (left column) and " << m_boundaryIndicator2
					<< " (right column):" << endl;
			std::map<double,
					std::pair<dealii::DoFHandler<2>::active_cell_iterator,
							size_t> >::iterator it1 = cellsAtBoundary1.begin();
			std::map<double,
					std::pair<dealii::DoFHandler<2>::cell_iterator, size_t> >::iterator it2 =
					cellsAtBoundary2.begin();
			while ((it1 != cellsAtBoundary1.end())
					or (it2 != cellsAtBoundary2.end())) {
				if (it1 != cellsAtBoundary1.end()) {
					info << it1->first << "  ";
					it1++;
				} else {
					info << "          ";
				}
				if (it2 != cellsAtBoundary2.end()) {
					info << it2->first;
					it2++;
				}
				info << endl;
			}
			throw PeriodicBoundaryNotPossible(
					"The discretizations of opposite periodic boundaries do not coincide. "
							"This version of the NATriuM solver does only work with equal discretizations. ",
					info);
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
	// Make iterators over active faces
	dealii::DoFHandler<3>::active_cell_iterator currentCell =
			doFHandler.begin_active();
	dealii::DoFHandler<3>::active_cell_iterator lastCell = doFHandler.end();
	// The key of these cells is a 2D point derived from the center point of each cell belonging to the boundary.
	// all the points in the respective boundaries have at least two coordinate components in common.  a dim-1 point
	//is derived from these common components and sorted in a way that a map can be made to connect the relative boundary cells.
	// the point_sorter class will enable us to order the key values in a specific fashion and the adjacent boundary cells
	//will be coupled to the same key values in the maps namely "cellsAtBoundary1" and "cellsAtBoundary2".
	// The second element of the value pair is the local face id of the face which belongs to the boundary.
	std::map<dealii::Point<2>,
			std::pair<dealii::DoFHandler<3>::active_cell_iterator, size_t>,
			BoundaryTools::point_sorter> cellsAtBoundary1;
	std::map<dealii::Point<2>,
			std::pair<dealii::DoFHandler<3>::cell_iterator, size_t>,
			BoundaryTools::point_sorter> cellsAtBoundary2;

	dealii::Point<2> key;

	// iterate over all active cells and sort them
	for (; currentCell != lastCell; ++currentCell) {
		if (currentCell->at_boundary()) {
			for (size_t i = 0; i < dealii::GeometryInfo<3>::faces_per_cell;
					i++) {
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator1) {
					//check if the face's normal vector is in X direction
					if ((dealii::GeometryInfo<3>::unit_normal_direction[i])
							== 0) {
						key[0] = (currentCell->face(i))->center()[1];
						key[1] = (currentCell->face(i))->center()[2];
					}
					// check if the face's normal vector is in Y direction
					else if (dealii::GeometryInfo<3>::unit_normal_direction[i]
							== 1) {
						key[0] = (currentCell->face(i))->center()[0];
						key[1] = (currentCell->face(i))->center()[2];
					} else {
						// if the normal vector is not in x or y direction then it must be in Z direction
						assert( (dealii::GeometryInfo<3>::unit_normal_direction[i])	== 2 );
						key[0] = (currentCell->face(i))->center()[0];
						key[1] = (currentCell->face(i))->center()[1];
					}

					cellsAtBoundary1.insert(
							std::make_pair(key,
									std::make_pair(currentCell,
											dealii::GeometryInfo<3>::opposite_face[i])));

				}
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator2) {
					if (dealii::GeometryInfo<3>::unit_normal_direction[i]
							== 0) {
						key[0] = (currentCell->face(i))->center()[1];
						key[1] = (currentCell->face(i))->center()[2];
					} else if (dealii::GeometryInfo<3>::unit_normal_direction[i]
							== 1) {
						key[0] = (currentCell->face(i))->center()[0];
						key[1] = (currentCell->face(i))->center()[2];
					} else {
						key[0] = (currentCell->face(i))->center()[0];
						key[1] = (currentCell->face(i))->center()[1];
					}
					cellsAtBoundary2.insert(
							std::make_pair(key,
									std::make_pair(currentCell,
											dealii::GeometryInfo<3>::opposite_face[i])));

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
	std::map<dealii::Point<2>,
			std::pair<dealii::DoFHandler<3>::active_cell_iterator, size_t> >::iterator atBoundary1 =
			cellsAtBoundary1.begin();
	std::map<dealii::Point<2>,
			std::pair<dealii::DoFHandler<3>::cell_iterator, size_t> >::iterator atBoundary2 =
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
		dealii::BlockDynamicSparsityPattern& cSparse, size_t n_blocks,
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
		dealii::DynamicSparsityPattern& block = cSparse.block(I, I);
		typename std::map<
				typename dealii::DoFHandler<dim>::active_cell_iterator,
				std::pair<typename dealii::DoFHandler<dim>::cell_iterator,
						size_t> >::const_iterator element = m_cells.begin();
		// for each cells belonging to the periodic boundary
		for (; element != m_cells.end(); element++) {
			element->first->get_dof_indices(doFIndicesAtCell1);
			element->second.first->get_dof_indices(doFIndicesAtCell2);
			// couple all dofs at boundary 1 with dofs at boundary 2
			constraints.add_entries_local_to_global(doFIndicesAtCell1,
					doFIndicesAtCell2, block, true);
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
		dealii::BlockDynamicSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const;
template void PeriodicBoundary<3>::addToSparsityPattern(
		dealii::BlockDynamicSparsityPattern& cSparse, size_t n_blocks,
		size_t n_dofs_per_block, size_t dofs_per_cell) const;

} /* namespace natrium */
