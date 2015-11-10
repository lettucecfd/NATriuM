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
#include "../utilities/Logging.h"

#include "BoundaryTools.h"

namespace natrium {

// The template Parameter has to be made explicit in order for the code to compile
template<size_t dim> PeriodicBoundary<dim>::PeriodicBoundary(
		size_t boundaryIndicator1, size_t boundaryIndicator2, size_t direction,
		shared_ptr<Mesh<dim> > triangulation) :
		m_boundaryIndicator1(boundaryIndicator1), m_boundaryIndicator2(
				boundaryIndicator2), m_direction(direction) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2) {
		throw PeriodicBoundaryNotPossible(
				"The boundary indicators defining the periodic boundaries must not be equal to each other.");
	}

	m_triangulation = triangulation;
	m_doFHandler = NULL;

	// collect periodic faces and extend halo layer of ghost cells across periodic boundaries
	// take care: collect_periodic_faces must not be called for a refined grid
	if (m_triangulation->n_levels() > 1) {
		throw PeriodicBoundaryNotPossible(
				"The periodic boundaries have to be set before the refinement.");
	}
	std::vector<
			dealii::GridTools::PeriodicFacePair<
					typename Mesh<dim>::cell_iterator> > matched_pairs;
	dealii::GridTools::collect_periodic_faces(*m_triangulation,
			m_boundaryIndicator1, m_boundaryIndicator2, m_direction,
			matched_pairs);
	// check if the boundary has neighbors
	if (matched_pairs.at(0).cell[0]->neighbor_index(
			matched_pairs.at(0).face_idx[0]) != -1) {
		throw PeriodicBoundaryNotPossible(
				"Cell at periodic boundary had neighbor. That is not allowed.");
	}
	// check if the boundary has already been created by setting user flags.
	// Somehow, this test does not work properly. The code crashes without
	// an error message in case you try to create a single Periodic Boundary
	// twice. I could not figure out why this does not work, but it should not
	// be so much of a problem, because normally the boundary collection will
	// already throw an error if you try to assign a boundary indicator to two
	// different objects in the boundary collection.
	for (size_t i = 0; i < matched_pairs.size(); i++) {
		size_t face_nr_1 = matched_pairs.at(i).face_idx[0];
		size_t face_nr_2 = matched_pairs.at(i).face_idx[1];
		if ((matched_pairs.at(i).cell[0]->face(face_nr_1)->user_flag_set())
				or (matched_pairs.at(i).cell[1]->face(face_nr_2)->user_flag_set())) {
			throw PeriodicBoundaryNotPossible(
					"You're trying to create a PeriodicBoundary, where another one "
							"has already been created. That does not work.");
		} else {
			matched_pairs.at(i).cell[0]->face(face_nr_1)->set_user_flag();
			matched_pairs.at(i).cell[1]->face(face_nr_2)->set_user_flag();
		}
	}

	LOG(DETAILED)
			<< "add periodicity to mesh... (if it succeeds, it will be displayed), in total "
			<< matched_pairs.size() << " matched pair(s) between boundaries "
			<< m_boundaryIndicator1 << " and  " << m_boundaryIndicator2 << endl;
	m_triangulation->add_periodicity(matched_pairs);
	LOG(DETAILED) << "...success" << endl;

} /* Constructor 2 */
template PeriodicBoundary<2>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2, size_t direction,
		shared_ptr<Mesh<2> > triangulation);
template PeriodicBoundary<3>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2, size_t direction,
		shared_ptr<Mesh<3> > triangulation);

template<size_t dim> bool PeriodicBoundary<dim>::isFaceInBoundary(
		const typename dealii::DoFHandler<dim>::active_cell_iterator &,
		size_t faceBoundaryIndicator) const {
	/*
	 // first condition: cell map has a key <cell>
	 if (m_cells.count(cell) == 0) {
	 return false;
	 }
	 // second condition: the face has the right boundary indicator
	 */
	// assert(m_cells.count(cell) != 0);
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
		const typename dealii::DoFHandler<2>::active_cell_iterator &,
		size_t faceBoundaryIndicator) const;
template bool PeriodicBoundary<3>::isFaceInBoundary(
		const typename dealii::DoFHandler<3>::active_cell_iterator &,
		size_t faceBoundaryIndicator) const;

template<size_t dim> PeriodicBoundary<dim>::~PeriodicBoundary() {
}
// The template Parameter has to be made explicit in order for the code to compile
template PeriodicBoundary<2>::~PeriodicBoundary();
template PeriodicBoundary<3>::~PeriodicBoundary();

template<size_t dim> void PeriodicBoundary<dim>::createCellMap(
		const dealii::DoFHandler<dim>& doFHandler) {
	m_doFHandler = &doFHandler;

	DealIIExtensions::make_periodicity_map_dg(doFHandler, m_boundaryIndicator1,
			m_boundaryIndicator2, m_direction, m_cells);

	checkCellMap();

} /* createMap */
template void PeriodicBoundary<2>::createCellMap(
		const dealii::DoFHandler<2>& doFHandler);
template void PeriodicBoundary<3>::createCellMap(
		const dealii::DoFHandler<3>& doFHandler);

template<size_t dim> void PeriodicBoundary<dim>::checkCellMap() {
	typename PeriodicCellMap<dim>::const_iterator it, end = m_cells.end();
	size_t face_nr_1, face_nr_2, boundary_id_1, boundary_id_2;
	for (it = m_cells.begin(); it != end; it++) {
		face_nr_1 = it->second.face_idx[0];
		face_nr_2 = it->second.face_idx[1];
		boundary_id_1 = it->second.cell[0]->face(face_nr_1)->boundary_id();
		boundary_id_2 = it->second.cell[1]->face(face_nr_2)->boundary_id();
		assert(
				(boundary_id_1 == m_boundaryIndicator1)
						or (boundary_id_1 == m_boundaryIndicator2));
		assert(
				(boundary_id_2 == m_boundaryIndicator1)
						or (boundary_id_2 == m_boundaryIndicator2));
		assert(not it->second.cell[0]->is_artificial());
		assert(not it->second.cell[1]->is_artificial());
		assert(it->second.cell[0]->active());
		assert(it->second.cell[1]->active());
		assert(
				(it->first == it->second.cell[0])
						or (it->first == it->second.cell[1]));

	}
}
template void PeriodicBoundary<2>::checkCellMap();
template void PeriodicBoundary<3>::checkCellMap();

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
	//initialize neighbor cell so that the dof handler is set properly,
	//in this way, m_cells can contain arbitrary TriaIterators and still return a DofAccessor
	//If this is left out we get an error in setting the neighbor cell
	neighborCell = m_doFHandler->begin_active();
	// return cell and face id
	FacePair<dim> face_pair = m_cells.at(cell);
	if (face_pair.cell[0] == cell) {
		neighborCell = face_pair.cell[1];
		return face_pair.face_idx[1];
	} else {
		assert(face_pair.cell[1] == cell);
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
		dealii::BlockDynamicSparsityPattern&, size_t, size_t, size_t) const {

	// THIS FUNCTION IS NOT USED!!! See DealIIExtensions module for details
	assert(false);
}
template void PeriodicBoundary<2>::addToSparsityPattern(
		dealii::BlockDynamicSparsityPattern&, size_t, size_t, size_t) const;
template void PeriodicBoundary<3>::addToSparsityPattern(
		dealii::BlockDynamicSparsityPattern&, size_t, size_t, size_t) const;

} /* namespace natrium */
