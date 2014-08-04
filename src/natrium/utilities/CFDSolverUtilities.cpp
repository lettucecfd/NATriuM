/**
 * @file CFDSolverUtilities.cpp
 * @short
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <utilities/CFDSolverUtilities.h>

#include "deal.II/base/geometry_info.h"

namespace natrium {

template<size_t dim>
double CFDSolverUtilities::getMinimumDoFDistanceGLL(
		const dealii::Triangulation<dim>& tria,
		const size_t orderOfFiniteElement) {
	assert(orderOfFiniteElement >= 2);
	// calculate minimal distance between vertices of the triangulation
	double min_vertex_distance = 100000000000.0;
	double distance = 0.0;
	for (typename dealii::Triangulation<dim>::active_cell_iterator cell =
			tria.begin_active(); cell != tria.end(); ++cell) {
		distance = cell->minimum_vertex_distance();
		if (distance < min_vertex_distance) {
			min_vertex_distance = distance;
		}
	}
	// calculate distance between closest quadrature nodes on a line
	dealii::QGaussLobatto<1> quadrature(orderOfFiniteElement);
	double min_dof_distance = 10000;
	for (size_t i = 0; i < orderOfFiniteElement; i++) {
		for (size_t j = i + 1; j < orderOfFiniteElement; j++) {
			double pointDist = quadrature.get_points().at(i).distance(
					quadrature.get_points().at(j));
			if (pointDist < min_dof_distance) {
				min_dof_distance = pointDist;
			}
		}
	}
	return min_vertex_distance * min_dof_distance;
}
template double CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
		const dealii::Triangulation<2>& tria, const size_t orderOfFiniteElement);
template double CFDSolverUtilities::getMinimumDoFDistanceGLL<3>(
		const dealii::Triangulation<3>& tria, const size_t orderOfFiniteElement);

} /* namespace natrium */
