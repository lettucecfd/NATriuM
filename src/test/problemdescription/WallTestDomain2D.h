/**
 * @file WallTestDomain2D.h
 * @short Description of a simple flow with walls (in square domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef WallTestDomain2D_H_
#define WallTestDomain2D_H_

#include "problemdescription/ProblemDescription.h"

#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_generator.h"

#include "problemdescription/MinLeeBoundary.h"
#include "utilities/BasicNames.h"

using dealii::Triangulation;

using namespace natrium;

/** @short Description of a simple Flow with wall boundaries (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of 4 elements.
 */
class WallTestDomain2D: public ProblemDescription<2> {
private:

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > makeGrid(size_t globalRefinementLevel) {

		//Creation of the principal domain
		shared_ptr<Triangulation<2> > unitSquare =
				make_shared<Triangulation<2> >();
		dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

		// Assign boundary indicators to the faces of the "parent cell"
		Triangulation<2>::active_cell_iterator cell =
				unitSquare->begin_active();
		cell->face(0)->set_all_boundary_indicators(0);  // left
		cell->face(1)->set_all_boundary_indicators(1);  // right
		cell->face(2)->set_all_boundary_indicators(2);  // top
		cell->face(3)->set_all_boundary_indicators(3);  // bottom

		// Refine grid to 2x2 = 4 cells
		unitSquare->refine_global(globalRefinementLevel);

		return unitSquare;
	}

	/**
	 * @short create boundaries
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries() {
		// make boundary description
		shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
				BoundaryCollection<2> >();
		boundaries->addBoundary(
				make_shared<MinLeeBoundary<2> >(0, numeric_vector(2)));
		boundaries->addBoundary(
				make_shared<MinLeeBoundary<2> >(1, numeric_vector(2)));
		boundaries->addBoundary(
				make_shared<MinLeeBoundary<2> >(2, numeric_vector(2)));
		numeric_vector topPlateVelocity(2);
		topPlateVelocity(0) = 0.01;
		boundaries->addBoundary(
				make_shared<MinLeeBoundary<2> >(3, topPlateVelocity));

		// Get the triangulation object (which belongs to the parent class).
		shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

		return boundaries;
	}

public:

	/// constructor
	WallTestDomain2D(size_t refineLevel):
		ProblemDescription<2>(makeGrid(refineLevel), 1.0, 1) {
			/// apply boundary values
			setBoundaries(makeBoundaries());
		}

	/// destructor
	virtual ~WallTestDomain2D(){};

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialDensities(distributed_vector& initialDensities,
			const vector<dealii::Point<2> >& supportPoints) const{
		for (size_t i = 0; i < initialDensities.size(); i++) {
			initialDensities(i) = 1.0;
		}
	}
	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<2> >& supportPoints) const{
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
			initialVelocities.at(0)(i) = 0.0;
			initialVelocities.at(1)(i) = 0.0;
		}
	}


};

#endif /* WallTestDomain2D_H_ */
