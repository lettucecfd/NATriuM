/**
 * @file CouetteFlow2D_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "benchmarks/CouetteFlow2D.h"

#include "boost/test/unit_test.hpp"

using std::endl;
using std::cout;

namespace natrium {


BOOST_AUTO_TEST_SUITE(CouetteFlow2D_test)

BOOST_AUTO_TEST_CASE(CouetteFlow2D_Construction_test) {
	cout << "CouetteFlow2D_Construction_test..." << endl;
	cout << "done" << endl;

} /* CouetteFlow2D_Construction_test */



BOOST_AUTO_TEST_CASE(CouetteFlow2D_Triangulation_test) {
	cout << "CouetteFlow2D_Triangulation_test..." << endl;
	cout << "done" << endl;

} /* CouetteFlow2D_Triangulation_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
