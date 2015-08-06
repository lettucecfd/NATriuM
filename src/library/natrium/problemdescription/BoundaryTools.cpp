/**
 * @file BoundaryTools.cpp
 * @short 
 * @date 14.11.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BoundaryTools.h"
#include "../utilities/Math.h"

bool natrium::BoundaryTools::checkParallelLines(
		const dealii::Point<2>& beginLine1, const dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2,
		std::string& errorMessage) {

	// check input
	double lengthLine1 = beginLine1.distance(endLine1);
	double lengthLine2 = beginLine2.distance(endLine2);

	// assert that points are not equal to one another
	if ((lengthLine1 == 0) or (lengthLine2 == 0)
			or (beginLine1.distance(beginLine2) == 0)
			or (endLine1.distance(endLine2) == 0)
			or (beginLine1.distance(endLine2) == 0)
			or (beginLine2.distance(endLine1) == 0)) {
		errorMessage.clear();
		errorMessage.append(
				"Two of the points defining a periodic boundary are equal. That is not allowed.");
		return false;
	}

	// assert that both lines have same length (up to 1%)
	if (abs(lengthLine1 - lengthLine2) / lengthLine1 > 0.01) {
		errorMessage.clear();
		errorMessage.append(
				"The two lines defining a periodic boundary must have the same length.");
		return false;
	}

	// assert that interfaces are parallel (anything else would need different handling)
	dealii::Tensor<1,2> differenceVector1 = endLine1 - beginLine1;
	dealii::Tensor<1,2> differenceVector2 = endLine2 - beginLine2;

	if (not Math::is_angle_small(differenceVector1, differenceVector2)) {
		// try to fix the problem by swapping begin and end
		dealii::Point<2> tmp = beginLine2;
		differenceVector2 = beginLine2 - endLine2;
		if (not Math::is_angle_small(differenceVector1, differenceVector2)) {
			// return
			errorMessage.clear();
			errorMessage.append(
					"The two lines defining a periodic boundary must be parallel to each other.");
			return false;
		}
		// re-orient line 2
		beginLine2 = endLine2;
		endLine2 = tmp;
	}
	return true;
}


bool natrium::BoundaryTools::getInterfacialLinesByBoundaryIndicator(
		size_t boundaryIndicator1, size_t boundaryIndicator2,
		shared_ptr<Mesh<2> > triangulation,
		dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2,
				std::string& errorMessage) {

	// Make iterators over active faces
	Mesh<2>::active_cell_iterator currentCell =
			triangulation->begin_active();
	Mesh<2>::active_cell_iterator lastCell =
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
		dealii::Tensor<1,2> line = endLine1 - beginLine1;
		dealii::Tensor<1,2> toPoint = element->second - beginLine1;
		if (not Math::is_angle_small(line, toPoint)) {
			std::stringstream s;
			s << "Not all points with boundary indicator "
					<< boundaryIndicator1 << " are on a line.";
			errorMessage = s.str();
			return false;
		}
	}
	for (element = ++pointsAtBoundary2.begin();
			element != --pointsAtBoundary2.end(); ++element) {
		// check if the vertex is really on line
		dealii::Tensor<1,2> line = endLine2 - beginLine2;
		dealii::Tensor<1,2> toPoint = element->second - beginLine2;
		if (not Math::is_angle_small(line, toPoint)) {
			std::stringstream s;
			s << "Not all points with boundary indicator "
					<< boundaryIndicator2 << " are on a line.";
			errorMessage = s.str();
			return false;
		}
	}

	return true;

}/* getInterfacialLinesByBoundaryIndicator */
