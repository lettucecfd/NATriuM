/**
 * @file CouetteFlow2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CouetteFlow2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

using dealii::GridGenerator;

namespace natrium {

CouetteFlow2D::CouetteFlow2D(double relaxationParameter,
		double topPlateVelocity) :
		ProblemDescription<2>(makeGrid(), relaxationParameter) {

	/// apply boundary values
	// setBoundaries(makeBoundaries(topPlateVelocity));

	/// set initial velocities to zero
	/// The numeric_vector is a dealii::Vector; The size-constructor applies the default value (0.0) to all components.
	numeric_vector initialVelocity(2);
	setConstantInitialVelocity(initialVelocity);

	/// set initial densities to 1.0
	setConstantInitialDensity(1.0);
}

CouetteFlow2D::~CouetteFlow2D() {
}

shared_ptr<Triangulation<2> > CouetteFlow2D::makeGrid() {

	//Creation of the principal domain
	shared_ptr<Triangulation<2> > unitSquare = make_shared<Triangulation<2> >();
	GridGenerator::hyper_cube(*unitSquare, 0, 1);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // bottom
	cell->face(1)->set_all_boundary_indicators(1);  //
	cell->face(2)->set_all_boundary_indicators(2);
	cell->face(3)->set_all_boundary_indicators(3);

	// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
	unitSquare->refine_global(3);

	return unitSquare;
}

shared_ptr<vector<BoundaryDescription<1> > > CouetteFlow2D::makeBoundaries(
		double topPlateVelocity) {

	// make boundary description
	shared_ptr<vector<BoundaryDescription<1> > > boundaries = make_shared<
			vector<BoundaryDescription<1> > >();

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}

} /* namespace natrium */
