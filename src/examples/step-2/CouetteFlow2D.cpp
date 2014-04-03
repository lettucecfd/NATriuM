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

#include "problemdescription/PeriodicBoundary.h"
#include "problemdescription/MinLeeBoundary.h"

namespace natrium {

CouetteFlow2D::CouetteFlow2D(double viscosity, double topPlateVelocity, size_t refinementLevel,
		double L) :
		ProblemDescription<2>(makeGrid(L, refinementLevel), viscosity, 1) {
	setCharacteristicLength(L);
	/// apply boundary values
	setBoundaries(makeBoundaries(topPlateVelocity));
}

CouetteFlow2D::~CouetteFlow2D() {
}

shared_ptr<Triangulation<2> > CouetteFlow2D::makeGrid(double L, size_t refinementLevel) {

	//Creation of the principal domain
	shared_ptr<Triangulation<2> > unitSquare = make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, L);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // bottom
	cell->face(3)->set_all_boundary_indicators(3);  // top

	// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
	unitSquare->refine_global(refinementLevel);

	return unitSquare;
}

shared_ptr<BoundaryCollection<2> > CouetteFlow2D::makeBoundaries(
		double topPlateVelocity) {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = topPlateVelocity;

	boundaries->addBoundary(
			make_shared<PeriodicBoundary<2> >(0, 1, getTriangulation()));
	boundaries->addBoundary(
			make_shared<MinLeeBoundary<2> >(2, zeroVelocity));
	boundaries->addBoundary(
				make_shared<MinLeeBoundary<2> >(3, constantVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}


double CouetteFlow2D::analyticVelocity1(const dealii::Point<2>& x,
		double t) const {
	// the analytic solution is given by an asymptotic series
	double U = getCharacteristicVelocity();
	double L = getCharacteristicLength();

	// the series converges veeeeeery slowly for t -> 0, thus assert t > epsilon
	assert (t > 0.1);

	double sum = U*x(1)/L;
	double lambda = 0.0;
	const double PI = atan(1)*4;
	double increment=1.0;
	for(size_t i = 1; fabs(increment) > 1e-12; i++){
		lambda = i*PI/L;
		increment = (i % 2 == 0 ? 1. : -1.) * 2*U/(lambda*L)*exp(-getViscosity()*lambda*lambda*t)*sin(lambda*x(1));
		// (i % 2 == 0 ? 1. : -1.) is a more efficient expression of (-1)^i =
		sum += increment;
		assert (i < 10000);
		//cout << i << " " << exp(-getViscosity()*lambda*lambda*t) << " " << sum << " " << (i % 2 == 0 ? 1. : -1.) * 2*U/(lambda*L)*exp(-getViscosity()*lambda*lambda*t)*sin(lambda*x(1)) << endl;
	}
	// assert convergence of the above asymptotic sum
	assert( fabs(increment) < 1e-9);
	return sum;
}

double CouetteFlow2D::analyticVelocity2(const dealii::Point<2>& x,
		double t) const {
	return 0.0;
}
} /* namespace natrium */
