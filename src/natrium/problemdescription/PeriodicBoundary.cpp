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

#include "../utilities/BasicNames.h"
#include "BoundaryTools.h"

namespace natrium {

/**
 * @short function to compare doubles as map keys;
 *        regards two doubles equal if they are in a epsilon(=1e-7)-range.
 */
class own_double_less: public std::binary_function<double, double, bool> {
public:
	own_double_less(double arg_ = 1e-7) :
			epsilon(arg_) {
	}
	bool operator()(const double &left, const double &right) const {
		// you can choose other way to make decision
		// (The original version is: return left < right;)
		return (abs(left - right) > epsilon) && (left < right);
	}
	double epsilon;
};

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
		getInterfacePositionsByBoundaryIndicator(boundaryIndicator1,
				boundaryIndicator2, triangulation, beginLine1, endLine1,
				beginLine2, endLine2);

		m_beginLine1 = beginLine1;
		m_beginLine2 = beginLine2;
		m_endLine1 = endLine1;
		m_endLine2 = endLine2;

		// check if positions of the interfaces are OK;
		// else throw PeriodicBoundaryNotPossible
		std::string errorMessage;
		bool isParallel = BoundaryTools::checkParallelLines(m_beginLine1,
				m_endLine1, m_beginLine2, m_endLine2, errorMessage);
		if (not isParallel)
			throw PeriodicBoundaryNotPossible(errorMessage);
	}
} /* Constructor 2 */
// The template Parameter has to be made explicit in order for the code to compile
template PeriodicBoundary<2>::PeriodicBoundary(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation);


template<size_t dim> bool PeriodicBoundary<dim>::isFaceInBoundary(const typename dealii::DoFHandler<dim>::active_cell_iterator & cell, size_t faceBoundaryIndicator) const {
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
template bool PeriodicBoundary<2>::isFaceInBoundary(const typename dealii::DoFHandler<2>::active_cell_iterator & cell, size_t faceBoundaryIndicator) const;
template bool PeriodicBoundary<3>::isFaceInBoundary(const typename dealii::DoFHandler<3>::active_cell_iterator & cell, size_t faceBoundaryIndicator) const;

template<size_t dim> PeriodicBoundary<dim>::~PeriodicBoundary() {
}
// The template Parameter has to be made explicit in order for the code to compile
template PeriodicBoundary<2>::~PeriodicBoundary();
template PeriodicBoundary<3>::~PeriodicBoundary();

template<> void PeriodicBoundary<2>::getInterfacePositionsByBoundaryIndicator(
		size_t boundaryIndicator1, size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation,
		dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2) {

	// Make iterators over active faces
	dealii::Triangulation<2>::active_cell_iterator currentCell =
			triangulation->begin_active();
	dealii::Triangulation<2>::active_cell_iterator lastCell =
			triangulation->end();

	// Make containers for all vertices at the boundary
	// maps are by default sorted by key;
	// The key of a point (x,y) is calculated 1000 * x + y
	// which leads to the fact that the first and last element of the map
	// define start and end point of the line.
	std::map<double, dealii::Point<2>, own_double_less> pointsAtBoundary1;
	std::map<double, dealii::Point<2>, own_double_less> pointsAtBoundary2;

	// iterate over all active faces and store coordinates of vertices in list
	for (; currentCell != lastCell; ++currentCell) {
		if (currentCell->at_boundary()) {
			for (size_t i = 0; i < dealii::GeometryInfo<2>::faces_per_cell;
					i++) {
				if (currentCell->face(i)->boundary_indicator()
						== boundaryIndicator1) {
					for (size_t j = 0;
							j < dealii::GeometryInfo<2>::vertices_per_face;
							j++) {
						double key = 1000. * currentCell->face(i)->vertex(j)[0]
								+ currentCell->face(i)->vertex(j)[1];
						pointsAtBoundary1.insert(
								std::make_pair(key,
										currentCell->face(i)->vertex(j)));
					}
				} else {
					if (currentCell->face(i)->boundary_indicator()
							== boundaryIndicator2) {
						for (size_t j = 0;
								j < dealii::GeometryInfo<2>::vertices_per_face;
								j++) {
							double key = 1000.
									* currentCell->face(i)->vertex(j)[0]
									+ currentCell->face(i)->vertex(j)[1];
							pointsAtBoundary2.insert(
									std::make_pair(key,
											currentCell->face(i)->vertex(j)));
						}
					}
				}
			}
		}
	}

	// get line edges
	beginLine1 = pointsAtBoundary1.begin()->second;
	endLine1 = (--pointsAtBoundary1.end())->second;
	beginLine2 = pointsAtBoundary2.begin()->second;
	endLine2 = (--pointsAtBoundary2.end())->second;

	// Make sure that vertices are really on the line
	// iterate over all active faces and store coordinates of vertices in list
	std::map<double, dealii::Point<2> >::iterator element;
	for (element = ++pointsAtBoundary1.begin();
			element != --pointsAtBoundary1.end(); ++element) {
		// check if the vertex is really on  line 1
		dealii::Point<2> line = endLine1 - beginLine1;
		dealii::Point<2> toPoint = element->second - beginLine1;
		if (not Math::is_angle_small(line, toPoint)) {
			std::stringstream errorMessage;
			errorMessage << "Not all points with boundary indicator "
					<< m_boundaryIndicator1 << " are on a line.";
			throw PeriodicBoundaryNotPossible(errorMessage.str());
		}
	}
	for (element = ++pointsAtBoundary2.begin();
			element != --pointsAtBoundary2.end(); ++element) {
		// check if the vertex is really on line
		dealii::Point<2> line = endLine2 - beginLine2;
		dealii::Point<2> toPoint = element->second - beginLine2;
		if (not Math::is_angle_small(line, toPoint)) {
			std::stringstream errorMessage;
			errorMessage << "Not all points with boundary indicator "
					<< m_boundaryIndicator2 << " are on a line.";
			throw PeriodicBoundaryNotPossible(errorMessage.str());
		}
	}

}/* getInterfacePositionsByBoundaryIndicator */

template<> void PeriodicBoundary<2>::createCellMap(
		const dealii::DoFHandler<2>& doFHandler) {

	// Make iterators over active faces
	dealii::DoFHandler<2>::active_cell_iterator currentCell =
			doFHandler.begin_active();
	dealii::DoFHandler<2>::active_cell_iterator lastCell = doFHandler.end();

	// The key of these cells is the distance from the begin point of the respective line.
	// The second element of the value pair is the local face id of the face which belongs to the boundary.
	std::map<double,
			std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t> > cellsAtBoundary1;
	std::map<double,
			std::pair<dealii::DoFHandler<2>::cell_iterator, size_t> > cellsAtBoundary2;

	// iterate over all active cells and sort them
	for (; currentCell != lastCell; ++currentCell) {
		if (currentCell->at_boundary()) {
			for (size_t i = 0; i < dealii::GeometryInfo<2>::faces_per_cell;
					i++) {
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator1) {
					double key = currentCell->center().distance(m_beginLine1);
					cellsAtBoundary1.insert(
							std::make_pair(key,
									std::make_pair(currentCell,
											dealii::GeometryInfo<2>::opposite_face[i])));
				}
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator2) {
					double key = currentCell->center().distance(m_beginLine2);
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
	std::map<double,
			std::pair<dealii::DoFHandler<2>::cell_iterator, size_t> >::iterator atBoundary2 =
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

template<> void PeriodicBoundary<2>::applyBoundaryValues(
		const shared_ptr<typename dealii::DoFHandler<2> > doFHandler,
		shared_ptr<dealii::ConstraintMatrix> constraintMatrix) const {

// check direction of the constraint
	size_t direction;
	if (std::fabs((m_beginLine1 - m_endLine1)[0]) < 1e-3) {
		direction = 0;
	} else if (std::fabs((m_beginLine1 - m_endLine1)[1]) < 1e-3) {
		direction = 1;
	} else {
		throw PeriodicBoundaryNotPossible(
				"The current implementation of periodic boundary conditions does only support boundaries parallel to the x- or y-axis.");
	}

	try {
// Note: For DoFHandler objects that are built on a parallel::distributed::Triangulation object
// parallel::distributed::Triangulation::add_periodicity has to be called before.
		dealii::DoFTools::make_periodicity_constraints(*doFHandler,
				m_boundaryIndicator1, m_boundaryIndicator2, direction,
				*constraintMatrix);
	} catch (std::exception& e) {
		throw PeriodicBoundaryNotPossible(
				"Error in dealii::DoFTools::make_periodicity_constraints");
	}

	/*
	 // Make iterators over active faces
	 dealii::DoFHandler<2>::active_cell_iterator currentCell =
	 doFHandler->begin_active();
	 size_t doFsPerCell = doFHandler->get_fe().dofs_per_cell;

	 // Sort dofs at boundary 1 with regard to their distance to beginLine1
	 std::map<double, size_t, own_double_less> doFsAtBoundary1;

	 // iterate over all active faces and store doFs of boundary 1 in list
	 for (; currentCell != doFHandler->end(); ++currentCell) {
	 if (currentCell->at_boundary()) {
	 feValues->reinit(currentCell);
	 for (size_t i = 0; i < dealii::GeometryInfo<2>::faces_per_cell; i++) {
	 if (currentCell->face(i)->boundary_indicator() == m_boundaryIndicator1) {

	 for (size_t j = 0;
	 j < feValues->dofs_per_face;
	 j++) {
	 durre
	 feValues->get_quadrature_points()
	 // TODO Iterate over non-vertex DOFs
	 double distance = m_beginLine1.distance(
	 currentCell->vertex(j));
	 doFsAtBoundary1.insert(
	 std::make_pair(distance,
	 currentCell->vertex_dof_index(j, 0)));
	 }
	 }
	 }
	 }
	 }

	 // iterate over all active faces and apply periodic boundary constraints to boundary 2
	 std::map<double, size_t>::iterator element;
	 for (; currentFace != doFHandler->end(); ++currentFace) {
	 if (currentFace->at_boundary()) {
	 if (currentFace->boundary_indicator() == m_boundaryIndicator2) {
	 for (size_t i = 0;
	 i < dealii::GeometryInfo<2>::vertices_per_face; i++) {
	 double distance = m_beginLine2.distance(
	 currentFace->vertex(i));

	 // get the respective dofs at boundary 1
	 // if key is already in map, element becomes the element with this key;
	 // and no insert takes place.
	 // if not: new element is inserted. Element points to new element
	 bool isKeyAlreadyInMap(false);
	 std::make_pair(element, isKeyAlreadyInMap) =
	 doFsAtBoundary1.insert(
	 std::make_pair(distance, 1e50));
	 if (isKeyAlreadyInMap) {
	 // add entry to constraint matrix which connects the doF with the respective
	 // doF at the other periodic interface
	 constraintMatrix->add_line(
	 currentFace->vertex_dof_index(i, 0));
	 constraintMatrix->add_entry(
	 currentFace->vertex_dof_index(i, 0),
	 element->second, 1.0);
	 } else {
	 // add entries to the constraint matrix which apply a linear mapping to the doF
	 // with the respective neighbors on the other periodic interface
	 double thisKey = element->first;
	 double smallerKey = (element--)->first;
	 double biggerKey = ((element++)++)->first;
	 constraintMatrix->add_line(
	 currentFace->vertex_dof_index(i, 0));
	 // element is now on this+1
	 constraintMatrix->add_entry(
	 currentFace->vertex_dof_index(i, 0),
	 element->second,
	 (thisKey - smallerKey)
	 / (biggerKey - smallerKey));
	 ((element--)--);
	 // element is now on this -1
	 constraintMatrix->add_entry(
	 currentFace->vertex_dof_index(i, 0),
	 element->second,
	 (biggerKey - thisKey)
	 / (biggerKey - smallerKey));
	 element++;
	 // erase element from list
	 assert(element->second > 1e49);
	 doFsAtBoundary1.erase(element);
	 }
	 }
	 }
	 }
	 }
	 */
} /* applyBoundaryValues */
template<> void PeriodicBoundary<3>::applyBoundaryValues(
		const shared_ptr<typename dealii::DoFHandler<3> > doFHandler,
		shared_ptr<dealii::ConstraintMatrix> constraintMatrix) const {
	//3D-Implementation
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
	return m_cells.at(cell).second;
}
// The template parameter has to be made explicit in order for the code to compile
template size_t PeriodicBoundary<2>::getOppositeCellAtPeriodicBoundary(
		const dealii::DoFHandler<2>::active_cell_iterator & cell,
		dealii::DoFHandler<2>::cell_iterator & neighborCell) const;
template size_t PeriodicBoundary<3>::getOppositeCellAtPeriodicBoundary(
		const dealii::DoFHandler<3>::active_cell_iterator & cell,
		dealii::DoFHandler<3>::cell_iterator & neighborCell) const;

} /* namespace natrium */
