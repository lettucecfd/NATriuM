/**
 * @file PeriodicTestDomain2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PeriodicTestDomain2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "problemdescription/PeriodicBoundary.h"

using namespace dealii::GridGenerator;

namespace natrium {

PeriodicTestDomain2D::PeriodicTestDomain2D(size_t globalRefinementLevel) :
		ProblemDescription<2>(makeGrid(globalRefinementLevel), 1.0) {

	/// apply boundary values
	setBoundaries(makeBoundaries());

}

PeriodicTestDomain2D::~PeriodicTestDomain2D() {
}

shared_ptr<Triangulation<2> > PeriodicTestDomain2D::makeGrid(size_t globalRefinementLevel) {

	//Creation of the principal domain
	shared_ptr<Triangulation<2> > unitSquare = make_shared<Triangulation<2> >();
	hyper_cube(*unitSquare, 0, 1);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // top
	cell->face(3)->set_all_boundary_indicators(3);  // bottom

	// Refine grid to 2x2 = 4 cells
	unitSquare->refine_global(globalRefinementLevel);

	return unitSquare;
}

shared_ptr<BoundaryCollection<2> > PeriodicTestDomain2D::makeBoundaries() {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	boundaries->addBoundary(make_shared<PeriodicBoundary<2> >(0,1,getTriangulation()));
	boundaries->addBoundary(make_shared<PeriodicBoundary<2> >(2,3,getTriangulation()));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}

} /* namespace natrium */
