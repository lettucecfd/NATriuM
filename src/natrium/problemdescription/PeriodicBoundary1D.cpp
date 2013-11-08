/**
 * @file PeriodicBoundary1D.cpp
 * @short Description of a periodic boundary on a line.
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PeriodicBoundary1D.h"

#include <iterator>

#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/dofs/dof_tools.h"

#include "../utilities/BasicNames.h"

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

/* DEPRECATED, because boundary indicators are needed for application of boundary values
 PeriodicBoundary1D::PeriodicBoundary1D(dealii::Point<2>& beginLine1,
 dealii::Point<2>& endLine1, dealii::Point<2>& beginLine2,
 dealii::Point<2>& endLine2,
 shared_ptr<dealii::Triangulation<2> > triangulation) {

 // check if positions of the interfaces are OK;
 // else throw PeriodicBoundaryNotPossible
 checkInterfacePositions(beginLine1, endLine1, beginLine2, endLine2);

 // check if lines really define a boundary of the triangulation
 checkInterfacesAtBoundary(beginLine1, endLine1, beginLine2, endLine2,
 triangulation);

 m_beginLine1 = beginLine1;
 m_beginLine2 = beginLine2;
 m_endLine1 = endLine1;
 m_endLine2 = endLine2;
 m_triangulation = triangulation;
 }*//* Constructor 1 */

PeriodicBoundary1D::PeriodicBoundary1D(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation) :
		m_boundaryIndicator1(boundaryIndicator1), m_boundaryIndicator2(
				boundaryIndicator2) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2) {
		throw PeriodicBoundaryNotPossible(
				"The boundary indicators defining the periodic boundaries must not be equal to each other.");
	}

	// create vertex points
	dealii::Point<2> beginLine1(0.0, 0.0);
	dealii::Point<2> endLine1(0.0, 0.0);
	dealii::Point<2> beginLine2(0.0, 0.0);
	dealii::Point<2> endLine2(0.0, 0.0);

	// calculate the positions of the vertex points(
	getInterfacePositionsByBoundaryIndicator(boundaryIndicator1,
			boundaryIndicator2, triangulation, beginLine1, endLine1, beginLine2,
			endLine2);

	m_beginLine1 = beginLine1;
	m_beginLine2 = beginLine2;
	m_endLine1 = endLine1;
	m_endLine2 = endLine2;
	m_triangulation = triangulation;

	// check if positions of the interfaces are OK;
	// else throw PeriodicBoundaryNotPossible
	checkInterfacePositions();

	// create the map which stores the opposite cells of each boundary cell, respectively
	createMap();

} /* Constructor 2 */

PeriodicBoundary1D::~PeriodicBoundary1D() {
}

void PeriodicBoundary1D::checkInterfacePositions() {

	// check input
	double lengthLine1 = m_beginLine1.distance(m_endLine1);
	double lengthLine2 = m_beginLine2.distance(m_endLine2);

	// assert that points are not equal to one another
	if ((lengthLine1 == 0) or (lengthLine2 == 0)
			or (m_beginLine1.distance(m_beginLine2) == 0)
			or (m_endLine1.distance(m_endLine2) == 0)
			or (m_beginLine1.distance(m_endLine2) == 0)
			or (m_beginLine2.distance(m_endLine1) == 0)) {
		throw PeriodicBoundaryNotPossible(
				"Two of the points defining a periodic boundary are equal. That is not allowed.");
	}

	// assert that both lines have same length (up to 1%)
	if (abs(lengthLine1 - lengthLine2) / lengthLine1 > 0.01) {
		throw PeriodicBoundaryNotPossible(
				"The two lines defining a periodic boundary must have the same length.");
	}

	// assert that interfaces are parallel (anything else would need different handling)
	dealii::Point<2> differenceVector1 = m_endLine1 - m_beginLine1;
	dealii::Point<2> differenceVector2 = m_endLine2 - m_beginLine2;
	if (not Math::is_angle_small(differenceVector1, differenceVector2)) {
		// try to fix the problem by swapping begin and end
		dealii::Point<2> tmp = m_beginLine2;
		m_beginLine2 = m_endLine2;
		m_endLine2 = tmp;
		differenceVector2 = m_endLine2 - m_beginLine2;
		if (not Math::is_angle_small(differenceVector1, differenceVector2)) {
			throw PeriodicBoundaryNotPossible(
					"The two lines defining a periodic boundary must be parallel to each other.");
		}
	}

} /* checkInterfacePositions */

void PeriodicBoundary1D::getInterfacePositionsByBoundaryIndicator(
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

void PeriodicBoundary1D::createMap() {

	// Make iterators over active faces
	dealii::Triangulation<2>::active_cell_iterator currentCell =
			m_triangulation->begin_active();
	dealii::Triangulation<2>::active_cell_iterator lastCell =
			m_triangulation->end();

	// The key of these cells is the distance from the begin point of the respective line.
	// The second element of the value pair is the local face id of the face which belongs to the boundary.
	std::map<double,
			std::pair<dealii::Triangulation<2>::active_cell_iterator, size_t> > cellsAtBoundary1;
	std::map<double,
			std::pair<dealii::Triangulation<2>::active_cell_iterator, size_t> > cellsAtBoundary2;

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
									std::make_pair(currentCell, dealii::GeometryInfo<2>::opposite_face[i])));
				}
				if (currentCell->face(i)->boundary_indicator()
						== m_boundaryIndicator2) {
					double key = currentCell->center().distance(m_beginLine2);
					cellsAtBoundary2.insert(
							std::make_pair(key,
									std::make_pair(currentCell, dealii::GeometryInfo<2>::opposite_face[i])));
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
	std::map<double,
			std::pair<dealii::Triangulation<2>::active_cell_iterator, size_t> >::iterator atBoundary1 =
			cellsAtBoundary1.begin();
	std::map<double,
			std::pair<dealii::Triangulation<2>::active_cell_iterator, size_t> >::iterator atBoundary2 =
			cellsAtBoundary2.begin();
	for (; atBoundary1 != --cellsAtBoundary1.end(); atBoundary1++) {
		dealii::CellId ID1 = atBoundary1->second.first->id();
		m_cells.insert(std::make_pair(ID1, atBoundary2->second));
		dealii::CellId ID2 = atBoundary2->second.first->id();
		m_cells.insert(std::make_pair(ID2, atBoundary1->second));
		atBoundary2++;
	}

} /* createMap */

void PeriodicBoundary1D::applyBoundaryValues(
		const shared_ptr<dealii::DoFHandler<2> > doFHandler,
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

size_t PeriodicBoundary1D::getOppositeCellAtPeriodicBoundary(
		dealii::CellId cellID,
		dealii::TriaIterator<dealii::CellAccessor<2> >& neighborCell) {

// assert that the given cell is at the boundary
	if (m_cells.count(cellID) == 0) {
		throw PeriodicBoundaryNotPossible(
				"The cellID does not belong to a cell at the boundary.");
	}

	neighborCell = m_cells[cellID].first;
	return m_cells[cellID].second;
}

} /* namespace natrium */
