/**
 * @file PoiseuilleFlow2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PoiseuilleFlow2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "../problemdescription/PeriodicBoundary.h"

#include "../utilities/Math.h"

namespace natrium {

PoiseuilleFlow2D::PoiseuilleFlow2D(double viscosity, size_t refinementLevel,
		double height, double length, bool is_periodic) :
		Benchmark<2>(makeGrid(refinementLevel, height, length, is_periodic),
				viscosity, height) {

	/// apply boundary values
	setBoundaries(makeBoundaries(is_periodic));

}

PoiseuilleFlow2D::~PoiseuilleFlow2D() {
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
shared_ptr<Triangulation<2> > PoiseuilleFlow2D::makeGrid(size_t refinementLevel,
		double height, double length, bool is_periodic) {
	//Creation of the principal domain
	shared_ptr<Triangulation<2> > square = make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_rectangle(*square, dealii::Point<2>(0, 0),
			dealii::Point<2>(1, 1), false);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell = square->begin_active();
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
shared_ptr<BoundaryCollection<2> > PoiseuilleFlow2D::makeBoundaries(bool is_periodic) {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	dealii::Vector<double> zeroVector(2);
	dealii::Vector<double> xVelocity(2);
	xVelocity(0) = 0.1 / sqrt(3);
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(0, zeroVector));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(1, zeroVector));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(2, zeroVector));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(3, xVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}
} /* namespace natrium */

void natrium::PoiseuilleFlow2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const {
}

void natrium::PoiseuilleFlow2D::getAnalyticVelocity(const dealii::Point<2>& x,
		double t, dealii::Point<2>& velocity) const {
}

double natrium::PoiseuilleFlow2D::getAnalyticDensity(const dealii::Point<2>& x,
		double t) const {
	return 1.0;
}
