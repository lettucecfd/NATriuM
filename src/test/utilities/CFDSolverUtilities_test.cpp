/**
 * @file CFDSolverUtilities_test.cpp
 * @short
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/CFDSolverUtilities.h"

#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_generator.h"

#include "boost/test/unit_test.hpp"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(CFDSolverUtilities_test)

BOOST_AUTO_TEST_CASE(CFDSolverUtilities_DoFDistance_test){
	pout << "CFDSolverUtilities_DoFDistance_test..." << endl;

	// make grid
	const double PI = 4*std::atan(1);
#ifdef WITH_TRILINOS_MPI
	Mesh<2> square(MPI_COMM_WORLD);
#else
	Mesh<2> square;
#endif
	dealii::GridGenerator::hyper_cube(square, 0.0, 2*PI);

	// test for different refinement levels
	square.refine_global(1);
	BOOST_CHECK_SMALL(2.0*PI/2.0 - CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(square, 1), 1e-10);
	square.refine_global(1);
	BOOST_CHECK_SMALL(2.0*PI/4.0 - CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(square, 1), 1e-10);

	// test for order 3 (dof in the middle of the cell)
	BOOST_CHECK_SMALL(2.0*PI/(4*2) - CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(square, 2), 1e-10);

	// test for order 4 (distances not regular any more)
	BOOST_CHECK_GT(2.0*PI/(4*3), CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(square, 3 ));

	pout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverUtilities_getIntegratorByID_test){
	pout << "CFDSolverUtilities_getIntegratorByID_test..." << endl;

	TimeIntegratorName t;
	DealIntegratorName d;
	string n;
	size_t id = 5;
	CFDSolverUtilities::get_integrator_by_id(id, t, d, n);
	BOOST_CHECK_EQUAL(t, OTHER);
	BOOST_CHECK_EQUAL(d, RK_THIRD_ORDER);
	BOOST_CHECK_EQUAL(n, "RK_THIRD_ORDER");

	pout << "done." << endl;
}


BOOST_AUTO_TEST_SUITE_END()

