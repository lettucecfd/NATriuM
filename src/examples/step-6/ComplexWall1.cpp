/**
 * @file ComplexWall1.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "ComplexWall1.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

#include "problemdescription/PeriodicBoundary.h"
#include "problemdescription/MinLeeBoundary.h"

#include "utilities/Logging.h"

namespace natrium {

ComplexWall1::ComplexWall1(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L) :
		ProblemDescription<2>(makeGrid(L, refinementLevel), viscosity, L), m_bottomVelocity(
				bottomVelocity) {
	setCharacteristicLength(L);

	/// apply boundary values
	setBoundaries (makeBoundaries(bottomVelocity));}

ComplexWall1::~ComplexWall1() {
}

shared_ptr<Triangulation<2> > ComplexWall1::makeGrid(double L,
		size_t refinementLevel) {

	//Creation of the principal domain
	shared_ptr<Triangulation<2> > unitSquare = make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, L);

	// Assign boundary indicators to the faces of the "parent cell"
	Triangulation<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // bottom
	cell->face(3)->set_all_boundary_indicators(3);  // top

	// refine grid
	unitSquare->refine_global(refinementLevel);

	// transform grid
		dealii::GridTools::transform(UnstructuredGridFunc(), *unitSquare);
		  std::ofstream out ("grid-2.eps");
		  dealii::GridOut grid_out;
		  grid_out.write_eps ( *unitSquare, out);
		return unitSquare;
}

shared_ptr<BoundaryCollection<2> > ComplexWall1::makeBoundaries(
		double bottomVelocity) {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = bottomVelocity;

	boundaries->addBoundary(
			make_shared<PeriodicBoundary<2> >(0, 1, getTriangulation()));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(2, constantVelocity));
	boundaries->addBoundary(
			make_shared<MinLeeBoundary<2> >(3, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

	return boundaries;
}

/**
 * @short set initial densities
 * @param[out] initialDensities vector of densities; to be filled
 * @param[in] supportPoints the coordinates associated with each degree of freedom
 */
void  ComplexWall1::applyInitialDensities(distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const{
	for (size_t i = 0; i < initialDensities.size(); i++){
		initialDensities(i) = 1.0;
	}
}

/**
 * @short set initial velocities
 * @param[out] initialVelocities vector of velocities; to be filled
 * @param[in] supportPoints the coordinates associated with each degree of freedom
 */
 void  ComplexWall1::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const{
	 assert ( initialVelocities.size() == 2);
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++){
			initialVelocities.at(0)(i) = 0.0;
			initialVelocities.at(1)(i) = 0.0;
		}
 }

} /* namespace natrium */
