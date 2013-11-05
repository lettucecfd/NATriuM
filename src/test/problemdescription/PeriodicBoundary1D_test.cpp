/**
 * @file PeriodicBoundary2D_test.cpp
 * @short 
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "problemdescription/PeriodicBoundary1D.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/base/point.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/lac/constraint_matrix.h"
#include "deal.II/fe/fe_q.h"

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(PeriodicBoundary1D_test)

/* DEPRECATED CONSTRUCTOR
 BOOST_AUTO_TEST_CASE(PeriodicBoundary1D_ConstructionByPoints_test) {

 cout << "PeriodicBoundary1D_ConstructionByPoints_test..." << endl;
 /////////////////
 // SANITY TEST //
 /////////////////
 /// Valid set of points and triangulation
 shared_ptr<dealii::Triangulation<2> > triangulation = make_shared<
 dealii::Triangulation<2> >();
 dealii::GridGenerator::hyper_cube(*triangulation, 0.0, 1.0);
 dealii::Point<2> beginLine1(0.0, 0.0);
 dealii::Point<2> endLine1(0.0, 1.0);
 dealii::Point<2> beginLine2(1.0, 0.0);
 dealii::Point<2> endLine2(1.0, 1.0);

 /// Test if construction of periodic boundary works
 BOOST_CHECK_NO_THROW(
 PeriodicBoundary1D(beginLine1, endLine1, beginLine2, endLine2,
 triangulation));

 /// Test if the same thing still works after refinement
 triangulation->refine_global(2);
 BOOST_CHECK_NO_THROW(
 PeriodicBoundary1D(beginLine1, endLine1, beginLine2, endLine2,
 triangulation));

 //////////////////
 // FAILURE TEST //
 //////////////////
 endLine2[1] += 0.1;
 BOOST_CHECK_THROW(
 PeriodicBoundary1D(beginLine1, endLine1, beginLine2, endLine2,
 triangulation), natrium::PeriodicBoundaryNotPossible);

 cout << "done." << endl;

 }
 */

BOOST_AUTO_TEST_CASE(PeriodicBoundary1D_ConstructionByBoundaryIndicator_test) {

	cout << "PeriodicBoundary1D_ConstructionByBoundaryIndicator_test..."
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
	BOOST_CHECK_NO_THROW(PeriodicBoundary1D(0, 1, triangulation));
	PeriodicBoundary1D myBoundary(0, 1, triangulation);
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
	BOOST_CHECK_NO_THROW(PeriodicBoundary1D(0, 1, triangulation));
	PeriodicBoundary1D myBoundary2(0, 1, triangulation);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine1()[0] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine1()[1] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine1()[0] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine1()[1] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine2()[0] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getBeginLine2()[1] - 0.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine2()[0] - 1.0, 1e-10);
	BOOST_CHECK_SMALL(myBoundary2.getEndLine2()[1] - 1.0, 1e-10);

	// Check if the same thing works for the top and bottom boundary
	BOOST_CHECK_NO_THROW(PeriodicBoundary1D(2, 3, triangulation));

	//////////////////
	// FAILURE TEST //
	//////////////////
	BOOST_CHECK_THROW(PeriodicBoundary1D(0, 0, triangulation),
			PeriodicBoundaryNotPossible);
	BOOST_CHECK_THROW(PeriodicBoundary1D(0, 2, triangulation),
			PeriodicBoundaryNotPossible);
	BOOST_CHECK_THROW(PeriodicBoundary1D(0, 3, triangulation),
			PeriodicBoundaryNotPossible);

	cout << "done." << endl;

} /* PeriodicBoundary1D_ConstructionByBoundaryIndicator_test */

BOOST_AUTO_TEST_CASE(PeriodicBoundary1D_ApplyBoundary_test) {

	cout << "PeriodicBoundary1D_ApplyBoundary_test..." << endl;

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
	PeriodicBoundary1D periodicLeftRight(0, 1, triangulation);
	PeriodicBoundary1D periodicTopBottom(2, 3, triangulation);

	// Apply boundary values
	periodicLeftRight.applyBoundaryValues(doFHandler, constraintMatrix);
	BOOST_CHECK(constraintMatrix->n_constraints() == std::pow(2,numberOfRefinementSteps)*(n_q_points)+1);
	periodicTopBottom.applyBoundaryValues(doFHandler, constraintMatrix);
	BOOST_CHECK(constraintMatrix->n_constraints() == 2*std::pow(2,numberOfRefinementSteps)*(n_q_points)+1);
	//constraintMatrix->print(cout);

	// Finalize constraint Matrix
	constraintMatrix->close();

	cout << "done." << endl;

} /* PeriodicBoundary1D_ApplyBoundary_test */

BOOST_AUTO_TEST_SUITE_END() /* PeriodicBoundary1D_test */

} /* namespace natrium */
