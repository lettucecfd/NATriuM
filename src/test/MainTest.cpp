/**
 * @file MainTest.cpp
 * @short 
 * @date 04.09.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main

#include "mpi.h"

#include "boost/test/unit_test.hpp"

#include "utilities/BasicNames.h"


using std::cout;
using std::endl;

struct MpiFixture {
	MpiFixture()   {
    	std::cout << "global setup\n";
#define WITH_TRILINOS
#ifdef WITH_TRILINOS
		// TRILINOS REQUIRES MPI TO BE INITIALIZED
		int argc = 0;
		char **argv;
		dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
#endif
    }
    ~MpiFixture()  { std::cout << "global teardown\n"; }
};

//____________________________________________________________________________//

BOOST_GLOBAL_FIXTURE( MpiFixture );


BOOST_AUTO_TEST_SUITE(Boost_test)


// Test if Boost unit test framework is running properly
BOOST_AUTO_TEST_CASE(Boost_test) {
	cout << "Boost_test..." << endl;
	BOOST_CHECK(1 == 1);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()
