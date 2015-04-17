/**
 * @file MainTest.cpp
 * @short 
 * @date 04.09.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main


#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/MPIGuard.h"


using std::cout;
using std::endl;


BOOST_AUTO_TEST_SUITE(Boost_test)


// Test if Boost unit test framework is running properly
BOOST_AUTO_TEST_CASE(Boost_test) {
	cout << "Boost_test..." << endl;

#ifdef WITH_TRILINOS
	natrium::MPIGuard::getInstance();
#endif

	BOOST_CHECK(1 == 1);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()
