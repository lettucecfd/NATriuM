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

#include "../problemdescription/PeriodicBoundary.h"

namespace natrium {

PeriodicTestDomain2D::PeriodicTestDomain2D(size_t globalRefinementLevel) :
		ProblemDescription<2>(makeGrid(), 1.0, 1) {

	/// apply boundary values
	setBoundaries(makeBoundaries());

	// Get the triangulation object (which belongs to the parent class).
	getMesh()->refine_global(globalRefinementLevel);


}

PeriodicTestDomain2D::~PeriodicTestDomain2D() {
}

boost::shared_ptr<Mesh<2> > PeriodicTestDomain2D::makeGrid() {

	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > unitSquare = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > unitSquare = boost::make_shared<Mesh<2> >();
#endif
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom

	// refinement has to be done after making the boundaries

	return unitSquare;
}

boost::shared_ptr<BoundaryCollection<2> > PeriodicTestDomain2D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(0,1,0,getMesh()));
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(2,3,1,getMesh()));



	return boundaries;
}

} /* namespace natrium */
