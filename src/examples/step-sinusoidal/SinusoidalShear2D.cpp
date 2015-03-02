/*
 * SinusoidalShear2D.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: kraemer
 */

#include "SinusoidalShear2D.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/grid_tools.h"

namespace natrium {

SinusoidalShear2D::SinusoidalShear2D(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L, double averageHeight,
		double amplitude) :
		ProblemDescription<2>(
				makeGrid(L, refinementLevel, averageHeight, amplitude),
				viscosity, L), m_bottomVelocity(bottomVelocity), m_height(averageHeight), m_ampl(amplitude) {
	setBoundaries(makeBoundaries(bottomVelocity));
}

SinusoidalShear2D::~SinusoidalShear2D() {
}

shared_ptr<Triangulation<2> > SinusoidalShear2D::makeGrid(double L,
		size_t refinementLevel, double averageHeight, double amplitude) {

	//Creation of the principal domain
	shared_ptr<Triangulation<2> > unitSquare = make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // bottom
	cell->face(3)->set_all_boundary_indicators(3);  // top

	// refine grid
	unitSquare->refine_global(refinementLevel);

	// transform grid
	dealii::GridTools::transform(UnstructuredGridFunc(averageHeight, amplitude, L),
			*unitSquare);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(*unitSquare, out);
	return unitSquare;
}

shared_ptr<BoundaryCollection<2> > SinusoidalShear2D::makeBoundaries(
		double bottomVelocity) {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = bottomVelocity;

	boundaries->addBoundary(
			make_shared<PeriodicBoundary<2> >(0, 1, getTriangulation()));
	boundaries->addBoundary(
			make_shared<MinLeeBoundary<2> >(2, constantVelocity));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(3, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}

void SinusoidalShear2D::applyInitialDensities(
		distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const {
	for (size_t i = 0; i < initialDensities.size(); i++) {
		initialDensities(i) = 1.0;
	}
}

void SinusoidalShear2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const {
	assert ( initialVelocities.size() == 2);
	for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
		double upper = m_height +  m_ampl * std::sin(8 * std::atan(1) * supportPoints.at(i)(0) );
		initialVelocities.at(0)(i) = m_bottomVelocity * pow( 1 - supportPoints.at(i)(1)/upper,2);
		initialVelocities.at(1)(i) = 0.0;
	}
}

} /* namespace natrium */
