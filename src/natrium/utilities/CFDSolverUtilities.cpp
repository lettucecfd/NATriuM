/**
 * @file CFDSolverUtilities.cpp
 * @short
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <utilities/CFDSolverUtilities.h>

#include <fstream>

#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/geometry_info.h"

namespace natrium {

template<size_t dim>
double CFDSolverUtilities::getMinimumDoFDistanceGLL(
		const dealii::Triangulation<dim>& tria,
		const size_t orderOfFiniteElement) {
	assert(orderOfFiniteElement >= 1);
	// calculate minimal distance between vertices of the triangulation
	double min_vertex_distance = CFDSolverUtilities::getMinimumVertexDistance<
			dim>(tria);

	// calculate distance between closest quadrature nodes on a line
	dealii::QGaussLobatto<1> quadrature(orderOfFiniteElement + 1);
	double min_dof_distance = 10000;
	for (size_t i = 0; i < orderOfFiniteElement + 1; i++) {
		for (size_t j = i + 1; j < orderOfFiniteElement + 1; j++) {
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
		const dealii::Triangulation<2>& tria,
		const size_t orderOfFiniteElement);
template double CFDSolverUtilities::getMinimumDoFDistanceGLL<3>(
		const dealii::Triangulation<3>& tria,
		const size_t orderOfFiniteElement);

template<size_t dim>
double CFDSolverUtilities::getMinimumVertexDistance(
		const dealii::Triangulation<dim>& tria) {
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
	if (min_vertex_distance < 1e-15) {
		throw CFDSolverUtilitiesException(
				"Minimum vertex distance < 1e-15; but the periodic boundary detection "
				"tests for cells distances < 1e-16 in own_double_less. That might cause "
				"serious problems.");
	}
	return min_vertex_distance;
}
template double CFDSolverUtilities::getMinimumVertexDistance<2>(
		const dealii::Triangulation<2>& tria);
template double CFDSolverUtilities::getMinimumVertexDistance<3>(
		const dealii::Triangulation<3>& tria);

template<size_t dim>
double CFDSolverUtilities::calculateTimestep(
		const dealii::Triangulation<dim>& tria,
		const size_t orderOfFiniteElement, const BoltzmannModel& boltzmannModel,
		double cFL) {
	assert(orderOfFiniteElement >= 1);
	double dx = CFDSolverUtilities::getMinimumVertexDistance<dim>(tria);
	double u = boltzmannModel.getMaxParticleVelocityMagnitude();
	// according to Hesthaven, dt ~ p^{-2}
	double dt = cFL * dx / (u * orderOfFiniteElement * orderOfFiniteElement);
	return dt;
}
template double CFDSolverUtilities::calculateTimestep<2>(
		const dealii::Triangulation<2>& tria, const size_t orderOfFiniteElement,
		const BoltzmannModel& boltzmannModel, double cFL);
template double CFDSolverUtilities::calculateTimestep<3>(
		const dealii::Triangulation<3>& tria, const size_t orderOfFiniteElement,
		const BoltzmannModel& boltzmannModel, double cFL);

template<int dim>
void CFDSolverUtilities::mesh_info(const dealii::Triangulation<dim> &tria,
		const std::string &filename) {
	std::cout << "Mesh info:" << std::endl << " dimension: " << dim << std::endl
			<< " no. of cells: " << tria.n_active_cells() << std::endl;
	{
		std::map<unsigned int, unsigned int> boundary_count;
		typename dealii::Triangulation<dim>::active_cell_iterator cell =
				tria.begin_active(), endc = tria.end();
		for (; cell != endc; ++cell) {
			for (unsigned int face = 0;
					face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
				if (cell->face(face)->at_boundary())
					boundary_count[cell->face(face)->boundary_indicator()]++;
			}
		}
		std::cout << " boundary indicators: ";
		for (std::map<unsigned int, unsigned int>::iterator it =
				boundary_count.begin(); it != boundary_count.end(); ++it) {
			std::cout << it->first << "(" << it->second << " times) ";
		}
		std::cout << std::endl;
	}
	std::ofstream out(filename.c_str());
	dealii::GridOut grid_out;
	grid_out.write_eps(tria, out);
	std::cout << " written to " << filename << std::endl << std::endl;
} /*mesh_info */
template void CFDSolverUtilities::mesh_info<2>(
		const dealii::Triangulation<2> &tria, const std::string &filename);
template void CFDSolverUtilities::mesh_info<3>(
		const dealii::Triangulation<3> &tria, const std::string &filename);

} /* namespace natrium */
