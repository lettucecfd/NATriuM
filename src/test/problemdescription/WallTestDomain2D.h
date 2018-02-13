/**
 * @file WallTestDomain2D.h
 * @short Description of a simple flow with walls (in square domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef WallTestDomain2D_H_
#define WallTestDomain2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_generator.h"

#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

using namespace natrium;

/** @short Description of a simple Flow with wall boundaries (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of 4 elements.
 */
class WallTestDomain2D: public ProblemDescription<2> {
private:

	size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid() {

		//Creation of the principal domain
		boost::shared_ptr<Mesh<2> > unitSquare = boost::make_shared<Mesh<2> >(
#ifdef WITH_TRILINOS_MPI
				MPI_COMM_WORLD
#endif
				);
		dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

		// Assign boundary indicators to the faces of the "parent cell"
		Mesh<2>::active_cell_iterator cell = unitSquare->begin_active();
		cell->face(0)->set_all_boundary_ids(0);  // left
		cell->face(1)->set_all_boundary_ids(1);  // right
		cell->face(2)->set_all_boundary_ids(2);  // top
		cell->face(3)->set_all_boundary_ids(3);  // bottom

		return unitSquare;
	}

	/**
	 * @short create boundaries
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries() {
		// make boundary description
		boost::shared_ptr<BoundaryCollection<2> > boundaries =
				boost::make_shared<BoundaryCollection<2> >();
		numeric_vector wall_velocity(2);
		wall_velocity(0) = 0.01;

		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
		/*boundaries->addBoundary(
				boost::make_shared<VelocityNeqBounceBack<2> >(0, wall_velocity));
		boundaries->addBoundary(
				boost::make_shared<VelocityNeqBounceBack<2> >(1, wall_velocity));
		*/
		boundaries->addBoundary(
				boost::make_shared<VelocityNeqBounceBack<2> >(2, wall_velocity));
		boundaries->addBoundary(
				boost::make_shared<VelocityNeqBounceBack<2> >(3, wall_velocity));

		return boundaries;
	}

public:

	/// constructor
	WallTestDomain2D(size_t refineLevel) :
			ProblemDescription<2>(makeGrid(), 0.2, 1), m_refinementLevel(
					refineLevel) {
		/// apply boundary values
		setBoundaries(makeBoundaries());
	}

	/// destructor
	virtual ~WallTestDomain2D() {
	}
	;


	virtual void refine(Mesh<2>& mesh) {
		// Refine grid to 2x2 = 4 cells
		mesh.refine_global(m_refinementLevel);
	}
	virtual void transform(Mesh<2>& ) {

	}
	virtual bool isCartesian(){
		return true;
	}

};

#endif /* WallTestDomain2D_H_ */
