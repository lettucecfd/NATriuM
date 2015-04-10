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

#include "../problemdescription/PeriodicBoundary.h"

#include "../utilities/Math.h"

namespace natrium {

TaylorGreenVortex2D::TaylorGreenVortex2D(double viscosity, size_t refinementLevel, double cs) :
		Benchmark<2>(makeGrid(refinementLevel), viscosity, 8*atan(1)), m_cs(cs){

	/// apply boundary values
	setBoundaries(makeBoundaries());

}

TaylorGreenVortex2D::~TaylorGreenVortex2D() {
}


void TaylorGreenVortex2D::getAnalyticVelocity(const dealii::Point<2>& x, double t, dealii::Point<2>& velocity) const {
	velocity(0) = sin(x(0))*cos(x(1))*exp(-2*getViscosity()*t);
	velocity(1) = -cos(x(0))*sin(x(1))*exp(-2*getViscosity()*t);
}


/**
 * @short get Analytic density at one point in space and time
 */
double TaylorGreenVortex2D::getAnalyticDensity(const dealii::Point<2>& x,
		double t) const {
	//double rho0 = 1;
	//double p = rho0/4.* (cos(2*x(0)) + cos(2*x(1))) * exp(-4 * getViscosity() * t);
	//return rho0 + p / (m_cs*m_cs) ;
	return 1.0;
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
shared_ptr<Triangulation<2> > TaylorGreenVortex2D::makeGrid(size_t refinementLevel) {
	//Creation of the principal domain
	shared_ptr<Triangulation<2> > square =
			make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*square, 0, 8*atan(1));

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
