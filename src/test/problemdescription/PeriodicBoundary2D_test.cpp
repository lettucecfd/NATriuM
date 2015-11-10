/**
 * @file PeriodicBoundary2D_test.cpp
 * @short 
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/problemdescription/PeriodicBoundary.h"

#include <iterator>

#include "boost/test/unit_test.hpp"

#include "deal.II/base/point.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/lac/constraint_matrix.h"
#include "deal.II/fe/fe_q.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"

#include "natrium/solver/SolverConfiguration.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

#include "TaylorGreenTest2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(PeriodicBoundary2D_test)

BOOST_AUTO_TEST_CASE(PeriodicBoundary2D_ConstructionByBoundaryIndicator_test) {

	pout << "PeriodicBoundary<2>_ConstructionByBoundaryIndicator_test..."
			<< endl;
	/////////////////
	// SANITY TEST //
	/////////////////
	/// Valid set of points and triangulation
	shared_ptr<Mesh<2> > triangulation = make_shared<Mesh<2> >(
#ifdef WITH_TRILINOS_MPI
			MPI_COMM_WORLD
#endif
			);
	dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);
	triangulation->begin_active()->face(0)->set_all_boundary_ids(0); //left
	triangulation->begin_active()->face(1)->set_all_boundary_ids(1); //right
	triangulation->begin_active()->face(2)->set_all_boundary_ids(2); //top
	triangulation->begin_active()->face(3)->set_all_boundary_ids(3); //bottom

	/// Test if construction works and vertices could be found automatically
	PeriodicBoundary<2>(0, 1, 0, triangulation);

	// Check if the same thing works for the top and bottom boundary
	BOOST_CHECK_NO_THROW(PeriodicBoundary<2>(2, 3, 1, triangulation));

	//////////////////
	// FAILURE TEST //
	//////////////////
	BOOST_CHECK_THROW(PeriodicBoundary<2>(0, 0, 0, triangulation),
			PeriodicBoundaryNotPossible);
	// Check that an error is thrown, when a boundary is created more than once
	// Ouch! This does not work, unfortunately -> see comment in PeriodicBoundary.cpp
	// BOOST_CHECK_THROW(PeriodicBoundary<2> myBoundary(0, 1, 0, triangulation), PeriodicBoundaryNotPossible);

	pout << "done." << endl;

} /* PeriodicBoundary<2>_ConstructionByBoundaryIndicator_test */

BOOST_AUTO_TEST_CASE(PeriodicBoundary2D_forDiscontinuousGalerkin_test) {

	pout << "PeriodicBoundary<2>_forDiscontinuousGalerkin_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	const size_t numberOfRefinementSteps = 2;

	/// create triangulation
	shared_ptr<Mesh<2> > triangulation = make_shared<Mesh<2> >(
#ifdef WITH_TRILINOS_MPI
			MPI_COMM_WORLD
#endif
			);
	dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);
	triangulation->begin_active()->face(0)->set_all_boundary_ids(0); //left
	triangulation->begin_active()->face(1)->set_all_boundary_ids(1); //right
	triangulation->begin_active()->face(2)->set_all_boundary_ids(2); //bottom
	triangulation->begin_active()->face(3)->set_all_boundary_ids(3); //top

	// make periodic boundaries object
	PeriodicBoundary<2> periodicLeftRight(0, 1, 0, triangulation);
	PeriodicBoundary<2> periodicTopBottom(2, 3, 1, triangulation);

	// refinement
	triangulation->refine_global(numberOfRefinementSteps);

	// distribute dofs
	shared_ptr<dealii::DoFHandler<2> > doFHandler = make_shared<
			dealii::DoFHandler<2> >(*triangulation);
	dealii::QGaussLobatto<1> qGaussLobatto(2);
	dealii::FE_DGQArbitraryNodes<2> fe(qGaussLobatto);
	doFHandler->distribute_dofs(fe);

	// create cell map
	periodicLeftRight.createCellMap(*doFHandler);
	periodicTopBottom.createCellMap(*doFHandler);

	// navigate to left upper corner
	dealii::DoFHandler<2>::active_cell_iterator leftUpperCorner =
			doFHandler->begin_active();
	bool cornerFound = false;
	for (; leftUpperCorner != doFHandler->end(); leftUpperCorner++) {
		for (size_t j = 0; j < 4; j++) {
			if ((leftUpperCorner->vertex(j)[0] < 0.00001)
					and (leftUpperCorner->vertex(j)[1] > 0.99999)) {
				cornerFound = true;
				break;
			}
		}
		if (cornerFound)
			break;
	}
	// make sure that the left upper corner was found
	BOOST_CHECK(cornerFound);
	BOOST_CHECK(leftUpperCorner->center()[0] == 0.125);
	BOOST_CHECK(leftUpperCorner->center()[1] == 0.875);
	BOOST_CHECK(leftUpperCorner->face(0)->boundary_id() == 0);
	BOOST_CHECK(leftUpperCorner->face(3)->boundary_id() == 3);

	// check the cell map defining the neighbors across boundaries
	const PeriodicCellMap<2>& cellMap = periodicLeftRight.getCellMap();
	/*std::map<dealii::DoFHandler<2>::active_cell_iterator,
	 std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t> >::const_iterator cell = cellMap.begin();
	 for (; cell != cellMap.end(); cell++){
	 pout << cell->first->center() <<  " -----" <<cell->second.second<<"----- "  << cell->second.first->center() << endl;

	 }*/
	BOOST_CHECK(cellMap.size() == 8);
	BOOST_CHECK(periodicTopBottom.getCellMap().size() == 8);

	// check function isFaceInBoundary (only left face can be in boundary)
	// for left upper cell
	BOOST_CHECK(
			periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(0)->boundary_id()));
	BOOST_CHECK(
			not periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(1)->boundary_id()));
	BOOST_CHECK(
			not periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(2)->boundary_id()));
	BOOST_CHECK(
			not periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(3)->boundary_id()));

	BOOST_CHECK(
			not periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(0)->boundary_id()));
	BOOST_CHECK(
			not periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(1)->boundary_id()));
	BOOST_CHECK(
			not periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(2)->boundary_id()));
	BOOST_CHECK(
			periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(3)->boundary_id()));

	// check if the opposite cells are really the opposite ones
	dealii::DoFHandler<2>::active_cell_iterator it, it2;
	size_t faceIndex = periodicLeftRight.getOppositeCellAtPeriodicBoundary(
			leftUpperCorner, it);
	size_t faceIndex2 = periodicTopBottom.getOppositeCellAtPeriodicBoundary(
			leftUpperCorner, it2);

	BOOST_CHECK(faceIndex == 1);
	BOOST_CHECK(faceIndex2 == 2);

	// Check if cells are correct
	BOOST_CHECK(it->face(1)->boundary_id() == 1);
	BOOST_CHECK(it->face(3)->boundary_id() == 3);
	BOOST_CHECK(it2->face(0)->boundary_id() == 0);
	BOOST_CHECK(it2->face(2)->boundary_id() == 2);

	//////////////////
	// FAILURE TEST //
	/////////////////

	// Not the same number of cells:
	leftUpperCorner->set_refine_flag();
	triangulation->execute_coarsening_and_refinement();
	//BOOST_CHECK_THROW(PeriodicBoundary<2>(0, 1, triangulation),
	//		PeriodicBoundaryNotPossible);

	doFHandler->clear();
	pout << "done." << endl;
} /*PeriodicBoundary<2>_forDisconitnuousGalerkin_test*/

BOOST_AUTO_TEST_SUITE_END() /* PeriodicBoundary<2>_test */

} /* namespace natrium */
