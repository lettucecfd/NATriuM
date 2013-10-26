/**
 * @file PeriodicBoundary1D.cpp
 * @short Description of a periodic boundary on a line.
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PeriodicBoundary1D.h"

#include <map>

#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

namespace natrium {

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
} /* Constructor 1 */

PeriodicBoundary1D::PeriodicBoundary1D(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation) {

	// check if boundary indcators are different
	if (boundaryIndicator1 == boundaryIndicator2){
		throw PeriodicBoundaryNotPossible("The boundary indicators defining the periodic boundaries must not be equal to each other.");
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

	/// SAME AS IN OTHER CONSTRUCTOR:
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
} /* Constructor 2 */

PeriodicBoundary1D::~PeriodicBoundary1D() {
}

void PeriodicBoundary1D::checkInterfacePositions(
		const dealii::Point<2>& beginLine1, const dealii::Point<2>& endLine1,
		const dealii::Point<2>& beginLine2,
		const dealii::Point<2>& endLine2) const {

	// check input
	double lengthLine1 = beginLine1.distance(endLine1);
	double lengthLine2 = beginLine2.distance(endLine2);

	// assert that points are not equal to one another
	if ((lengthLine1 == 0) or (lengthLine2 == 0)
			or (beginLine1.distance(beginLine2) == 0)
			or (endLine1.distance(endLine2) == 0)
			or (beginLine1.distance(endLine2) == 0)
			or (beginLine2.distance(endLine1) == 0)) {
		throw PeriodicBoundaryNotPossible(
				"Two of the points defining a periodic boundary are equal. That is not allowed.");
	}

	// assert that both lines have same length (up to 1%)
	if (abs(lengthLine1 - lengthLine2) / lengthLine1 > 0.01) {
		throw PeriodicBoundaryNotPossible(
				"The two lines defining a periodic boundary must have the same length.");
	}

	// assert that interfaces are parallel (anything else would need different handling)
	dealii::Point<2> differenceVector1 = endLine1 - beginLine1;
	dealii::Point<2> differenceVector2 = endLine2 - beginLine2;
	if (abs(differenceVector1 * differenceVector2)
			/ (differenceVector1.norm() * differenceVector2.norm()) < 0.99) {
		throw PeriodicBoundaryNotPossible(
				"The two lines defining a periodic boundary must be parallel to each other.");
	}

} /* checkInterfacePositions */

void PeriodicBoundary1D::checkInterfacesAtBoundary(
		const dealii::Point<2>& beginLine1, const dealii::Point<2>& endLine1,
		const dealii::Point<2>& beginLine2, const dealii::Point<2>& endLine2,
		shared_ptr<dealii::Triangulation<2> > triangulation) {
} /* checkInterfacesAtBoundary */

void PeriodicBoundary1D::getInterfacePositionsByBoundaryIndicator(
		size_t boundaryIndicator1, size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation,
		dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2) {

	// Make iterators over active faces
	dealii::Triangulation<2>::active_face_iterator currentFace =
			triangulation->begin_active_face();
	dealii::Triangulation<2>::active_face_iterator lastFace =
			triangulation->end_face();

	// Make containers for all vertices at the boundary
	// maps are by default sorted by key;
	// The key of a point (x,y) is calculated 1000 * x + y
	// which leads to the fact that the first and last element of the map
	// define start and end point of the line.
	std::map<double, dealii::Point<2> > pointsAtBoundary1;
	std::map<double, dealii::Point<2> > pointsAtBoundary2;

	// iterate over all active faces and store coordinates of vertices in list
	for (; currentFace != lastFace; ++currentFace) {
		if (currentFace->at_boundary()) {
			if (currentFace->boundary_indicator() == boundaryIndicator1) {
				for (size_t i = 0;
						i < dealii::GeometryInfo<2>::vertices_per_face; i++) {
					double key = 1000. * currentFace->vertex(i)[0]
							+ currentFace->vertex(i)[1];
					/*cout << "vertex to 1: " << currentFace->vertex(i)[0] << " "
							<< currentFace->vertex(i)[1] << "; key: " << key
							<< endl;*/
					pointsAtBoundary1.insert(
							std::make_pair(key, currentFace->vertex(i)));
				}
			} else {
				if (currentFace->boundary_indicator() == boundaryIndicator2) {
					for (size_t i = 0;
							i < dealii::GeometryInfo<2>::vertices_per_face;
							i++) {
						double key = 1000. * currentFace->vertex(i)[0]
								+ currentFace->vertex(i)[1];
						/*cout << "vertex to 2: " << currentFace->vertex(i)[0] << " "
								<< currentFace->vertex(i)[1] << "; key: " << key
								<< endl;*/
						pointsAtBoundary2.insert(
								std::make_pair(key, currentFace->vertex(i)));
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

} /* getInterfacePositionsByBoundaryIndicator */

void PeriodicBoundary1D::applyBoundaryValues(
		shared_ptr<dealii::DoFHandler<2> > doFHandler) {

// assert that points are at faces of the mesh

} /* applyBoundaryValues */

} /* namespace natrium */
