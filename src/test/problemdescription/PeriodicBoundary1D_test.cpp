/**
 * @file PeriodicBoundary2D_test.cpp
 * @short 
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "problemdescription/PeriodicBoundary.h"

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

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(PeriodicBoundary2D_test)


BOOST_AUTO_TEST_CASE(PeriodicBoundary2D_ConstructionByBoundaryIndicator_test) {

	cout << "PeriodicBoundary<2>_ConstructionByBoundaryIndicator_test..."
			<< endl;
	/////////////////
	// SANITY TEST //
	/////////////////
	/// Valid set of points and triangulation
	shared_ptr<dealii::Triangulation<2> > triangulation = make_shared<
			dealii::Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);
	triangulation->begin_active(0)->face(0)->set_boundary_indicator(0); //left
	triangulation->begin_active(0)->face(1)->set_boundary_indicator(1); //right
	triangulation->begin_active(0)->face(2)->set_boundary_indicator(2); //top
	triangulation->begin_active(0)->face(3)->set_boundary_indicator(3); //bottom

	/// Test if construction works and vertices could be found automatically
	BOOST_CHECK_NO_THROW(PeriodicBoundary<2>(0, 1, triangulation));
	PeriodicBoundary<2> myBoundary(0, 1, triangulation);
	BOOST_CHECK_SMALL(myBoundary.getBeginLine1()[0] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getBeginLine1()[1] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getEndLine1()[0] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getEndLine1()[1] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getBeginLine2()[0] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getBeginLine2()[1] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getEndLine2()[0] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary.getEndLine2()[1] - 1.0, 1e-10);

	// Check if same thing still works after refinement
	triangulation->refine_global(2);
	BOOST_CHECK_NO_THROW(PeriodicBoundary<2>(0, 1, triangulation));
	PeriodicBoundary<2> myBoundary2(0, 1, triangulation);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine1()[0] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine1()[1] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine1()[0] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine1()[1] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine2()[0] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine2()[1] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine2()[0] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine2()[1] - 1.0, 1e-10);

	// Check if the same thing works for the top and bottom boundary
	BOOST_CHECK_NO_THROW(PeriodicBoundary<2>(2, 3, triangulation));

	//////////////////
	// FAILURE TEST //
	//////////////////
	BOOST_CHECK_THROW(PeriodicBoundary<2>(0, 0, triangulation),
			PeriodicBoundaryNotPossible);
	BOOST_CHECK_THROW(PeriodicBoundary<2>(0, 2, triangulation),
			PeriodicBoundaryNotPossible);
	BOOST_CHECK_THROW(PeriodicBoundary<2>(0, 3, triangulation),
			PeriodicBoundaryNotPossible);

	cout << "done." << endl;

} /* PeriodicBoundary<2>_ConstructionByBoundaryIndicator_test */

BOOST_AUTO_TEST_CASE(PeriodicBoundary2D_ApplyBoundary_test) {

	cout << "PeriodicBoundary<2>_ApplyBoundary_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	const size_t numberOfRefinementSteps = 3;
	const size_t n_q_points = 2;

	/// create triangulation
	shared_ptr<dealii::Triangulation<2> > triangulation = make_shared<
			dealii::Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);
	triangulation->begin_active(0)->face(0)->set_boundary_indicator(0); //left
	triangulation->begin_active(0)->face(1)->set_boundary_indicator(1); //right
	triangulation->begin_active(0)->face(2)->set_boundary_indicator(2); //top
	triangulation->begin_active(0)->face(3)->set_boundary_indicator(3); //bottom
	triangulation->refine_global(numberOfRefinementSteps);

	// create finite element, dof handler and constraint matrix
	shared_ptr<dealii::FE_Q<2> > fe = make_shared<dealii::FE_Q<2> >(n_q_points);
	shared_ptr<dealii::DoFHandler<2> > doFHandler = make_shared<
			dealii::DoFHandler<2> >(*triangulation);
	shared_ptr<dealii::ConstraintMatrix> constraintMatrix = make_shared<
			dealii::ConstraintMatrix>();
	constraintMatrix->clear();

	// distribute degrees of freedom
	doFHandler->distribute_dofs(*fe);

	// make periodic boundaries object
	PeriodicBoundary<2> periodicLeftRight(0, 1, triangulation);
	PeriodicBoundary<2> periodicTopBottom(2, 3, triangulation);

	// Apply boundary values
	periodicLeftRight.applyBoundaryValues(doFHandler, constraintMatrix);
	BOOST_CHECK(
			constraintMatrix->n_constraints()
					== std::pow(2, numberOfRefinementSteps) * (n_q_points) + 1);
	periodicTopBottom.applyBoundaryValues(doFHandler, constraintMatrix);
	BOOST_CHECK(
			constraintMatrix->n_constraints()
					== 2 * std::pow(2, numberOfRefinementSteps) * (n_q_points)
							+ 1);
	//constraintMatrix->print(cout);

	// Finalize constraint Matrix
	constraintMatrix->close();

	cout << "done." << endl;

} /* PeriodicBoundary<2>_ApplyBoundary_test */

BOOST_AUTO_TEST_CASE(PeriodicBoundary2D_forDiscontinuousGalerkin_test) {

	cout << "PeriodicBoundary<2>_forDiscontinuousGalerkin_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	const size_t numberOfRefinementSteps = 2;

	/// create triangulation
	shared_ptr<dealii::Triangulation<2> > triangulation = make_shared<
			dealii::Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);
	triangulation->begin_active(0)->face(0)->set_boundary_indicator(0); //left
	triangulation->begin_active(0)->face(1)->set_boundary_indicator(1); //right
	triangulation->begin_active(0)->face(2)->set_boundary_indicator(2); //bottom
	triangulation->begin_active(0)->face(3)->set_boundary_indicator(3); //top
	triangulation->refine_global(numberOfRefinementSteps);

	// make periodic boundaries object
	PeriodicBoundary<2> periodicLeftRight(0, 1, triangulation);
	PeriodicBoundary<2> periodicTopBottom(2, 3, triangulation);

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
	BOOST_CHECK(leftUpperCorner->face(0)->boundary_indicator() == 0);
	BOOST_CHECK(leftUpperCorner->face(3)->boundary_indicator() == 3);

	// check the cell map defining the neighbors across boundaries
	const std::map<dealii::DoFHandler<2>::active_cell_iterator,
				std::pair<dealii::DoFHandler<2>::cell_iterator, size_t> >& cellMap = periodicLeftRight.getCellMap();
	/*std::map<dealii::DoFHandler<2>::active_cell_iterator,
					std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t> >::const_iterator cell = cellMap.begin();
	for (; cell != cellMap.end(); cell++){
		cout << cell->first->center() <<  " -----" <<cell->second.second<<"----- "  << cell->second.first->center() << endl;

	}*/
	BOOST_CHECK(cellMap.size()==8);
	BOOST_CHECK(periodicTopBottom.getCellMap().size() == 8);

	// check function isFaceInBoundary (only left face can be in boundary)
	// for left upper cell
	BOOST_CHECK(
			periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(0)->boundary_indicator()));
	BOOST_CHECK(
			not periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(1)->boundary_indicator()));
	BOOST_CHECK(
			not periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(2)->boundary_indicator()));
	BOOST_CHECK(
			not periodicLeftRight.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(3)->boundary_indicator()));

	BOOST_CHECK(
			not periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(0)->boundary_indicator()));
	BOOST_CHECK(
			not periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(1)->boundary_indicator()));
	BOOST_CHECK(
			not periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(2)->boundary_indicator()));
	BOOST_CHECK(
			periodicTopBottom.isFaceInBoundary(leftUpperCorner,
					leftUpperCorner->face(3)->boundary_indicator()));

	// check if the opposite cells are really the opposite ones
	dealii::DoFHandler<2>::active_cell_iterator it, it2;
	size_t faceIndex = periodicLeftRight.getOppositeCellAtPeriodicBoundary(
			leftUpperCorner, it);
	size_t faceIndex2 = periodicTopBottom.getOppositeCellAtPeriodicBoundary(
			leftUpperCorner, it2);

	BOOST_CHECK(faceIndex == 0);
	BOOST_CHECK(faceIndex2 == 3);

	// Check if cells are correct
	BOOST_CHECK(it->face(1)->boundary_indicator() == 1);
	BOOST_CHECK(it->face(3)->boundary_indicator() == 3);
	BOOST_CHECK(it2->face(0)->boundary_indicator() == 0);
	BOOST_CHECK(it2->face(2)->boundary_indicator() == 2);

	//////////////////
	// FAILURE TEST //
	/////////////////

	// Not the same number of cells:
	leftUpperCorner->set_refine_flag();
	triangulation->execute_coarsening_and_refinement();
	//BOOST_CHECK_THROW(PeriodicBoundary<2>(0, 1, triangulation),
	//		PeriodicBoundaryNotPossible);

	doFHandler->clear();
	cout << "done." << endl;
} /*PeriodicBoundary<2>_forDisconitnuousGalerkin_test*/

BOOST_AUTO_TEST_SUITE_END() /* PeriodicBoundary<2>_test */

} /* namespace natrium */
