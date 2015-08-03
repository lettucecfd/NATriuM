/**
 * @file LidDrivenCavity2D.cpp
 * @short Lid-driven cavity with three static walls and one moving wall
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "LidDrivenCavity2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "natrium/problemdescription/PeriodicBoundary.h"

#include "natrium/utilities/Math.h"

namespace natrium {

LidDrivenCavity2D::LidDrivenCavity2D(double velocity, double viscosity,
		size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(refinementLevel), viscosity, 1.0), topPlateVelocity(
				velocity) {

	/// apply boundary values
	setBoundaries(makeBoundaries());

}

LidDrivenCavity2D::~LidDrivenCavity2D() {
}

void LidDrivenCavity2D::applyInitialDensities(
		distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const {
	for (size_t i = 0; i < initialDensities.size(); i++) {
		initialDensities(i) = 1.0;
	}
}

void LidDrivenCavity2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const {
	assert(initialVelocities.at(0).size() == initialVelocities.at(1).size());
	for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
		initialVelocities.at(0)(i) = 0.0;
		initialVelocities.at(1)(i) = 0.0;
	}
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
shared_ptr<Mesh<2> > LidDrivenCavity2D::makeGrid(
		size_t refinementLevel) {
	//Creation of the principal domain
	shared_ptr<Mesh<2> > square = make_shared<Mesh<2> >();
	dealii::GridGenerator::hyper_rectangle(*square, dealii::Point<2>(0, 0),
			dealii::Point<2>(1, 1), false);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = square->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // bottom
	cell->face(3)->set_all_boundary_indicators(3);  // top

	// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
	square->refine_global(refinementLevel);

	return square;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
shared_ptr<BoundaryCollection<2> > LidDrivenCavity2D::makeBoundaries() {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	dealii::Vector<double> zeroVector(2);
	dealii::Vector<double> xVelocity(2);
	xVelocity(0) = topPlateVelocity;
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(0, zeroVector));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(1, zeroVector));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(2, zeroVector));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(3, xVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
