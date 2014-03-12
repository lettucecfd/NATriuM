/**
 * @file TaylorGreenVortex2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TaylorGreenVortex2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "problemdescription/PeriodicBoundary.h"

#include "utilities/Math.h"

namespace natrium {

TaylorGreenVortex2D::TaylorGreenVortex2D(double viscosity, size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(refinementLevel), viscosity, 1) {

	/// apply boundary values
	setBoundaries(makeBoundaries());

}

TaylorGreenVortex2D::~TaylorGreenVortex2D() {
}


void TaylorGreenVortex2D::applyInitialDensities(
		distributed_vector& initialDensities,
		vector<dealii::Point<2> >& supportPoints) const {
	for (size_t i = 0; i < initialDensities.size(); i++) {
		initialDensities(i) = 1.0;
	}
}

void TaylorGreenVortex2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		vector<dealii::Point<2> >& supportPoints) const {
	assert(
			initialVelocities.at(0).size()
					== initialVelocities.at(1).size());
	for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
		initialVelocities.at(0)(i) = analyticVelocity1(supportPoints.at(i), 0);
		initialVelocities.at(1)(i) = analyticVelocity2(supportPoints.at(i), 0);;
	}
}

double TaylorGreenVortex2D::analyticVelocity1(const dealii::Point<2>& x,
		double t) const {
	return sin(x(0))*cos(x(1))*exp(-2*getViscosity()*t);
}

double TaylorGreenVortex2D::analyticVelocity2(const dealii::Point<2>& x,
		double t) const {
	return -cos(x(0))*sin(x(1))*exp(-2*getViscosity()*t);
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
shared_ptr<Triangulation<2> > TaylorGreenVortex2D::makeGrid(size_t refinementLevel) {
	//Creation of the principal domain
	shared_ptr<Triangulation<2> > square =
			make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*square, 0, 2*Math::PI);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell =
			square->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // top
	cell->face(3)->set_all_boundary_indicators(3);  // bottom

	// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
	square->refine_global(refinementLevel);

	return square;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
shared_ptr<BoundaryCollection<2> > TaylorGreenVortex2D::makeBoundaries() {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	boundaries->addBoundary(
			make_shared<PeriodicBoundary<2> >(0, 1, getTriangulation()));
	boundaries->addBoundary(
			make_shared<PeriodicBoundary<2> >(2, 3, getTriangulation()));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}
} /* namespace natrium */
