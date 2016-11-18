/**
 * @file MainTest.cpp
 * @short 
 * @date 04.09.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main


#include <boost/test/included/unit_test.hpp>

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/MPIGuard.h"

using namespace natrium;

/* MPI has to be initialized before MPI_COMM_WORLD is duplicated
 * That happens quite often in the code, so we have to make sure
 * that MPIInitFinalize has been called before the tests start.
 * This is done by a global fixture:
 */
struct InitMPI {
	InitMPI()   {
    	std::cout << "Global fixture: MPI Init\n" << endl;
    	MPIGuard::getInstance();
    }
    ~InitMPI()  {  }
};

BOOST_AUTO_TEST_SUITE(Boost_test)


BOOST_GLOBAL_FIXTURE( InitMPI );


// Test if Boost unit test framework is running properly
BOOST_AUTO_TEST_CASE(Boost_test) {
	pout << "Boost_test..." << endl;

	BOOST_CHECK(1 == 1);
	pout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

