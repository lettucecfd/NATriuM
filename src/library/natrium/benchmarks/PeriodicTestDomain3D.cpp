/**
 * @file PeriodicTestDomain3D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PeriodicTestDomain3D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "../problemdescription/PeriodicBoundary.h"

namespace natrium {

PeriodicTestDomain3D::PeriodicTestDomain3D(size_t globalRefinementLevel) :
		ProblemDescription<3>(makeGrid(globalRefinementLevel), 1.0, 1) {

	/// apply boundary values
	setBoundaries(makeBoundaries());

}

PeriodicTestDomain3D::~PeriodicTestDomain3D() {
}

shared_ptr<Mesh<3> > PeriodicTestDomain3D::makeGrid(size_t globalRefinementLevel) {

	//Creation of the principal domain
	shared_ptr<Mesh<3> > unitSquare = make_shared<Mesh<3> >();
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<3>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // fron
	cell->face(3)->set_all_boundary_indicators(3);  // back
	cell->face(4)->set_all_boundary_indicators(4);  // bottom
	cell->face(5)->set_all_boundary_indicators(5);  // top

	// Refine grid to 2x2 = 4 cells
	unitSquare->refine_global(globalRefinementLevel);

	return unitSquare;
}

shared_ptr<BoundaryCollection<3> > PeriodicTestDomain3D::makeBoundaries() {

	// make boundary description
	shared_ptr<BoundaryCollection<3> > boundaries = make_shared<
			BoundaryCollection<3> >();
	boundaries->addBoundary(make_shared<PeriodicBoundary<3> >(0,1,getMesh()));
	boundaries->addBoundary(make_shared<PeriodicBoundary<3> >(2,3,getMesh()));
	boundaries->addBoundary(make_shared<PeriodicBoundary<3> >(4,5,getMesh()));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Mesh<3> > tria_pointer = getMesh();

	return boundaries;
}

} /* namespace natrium */
