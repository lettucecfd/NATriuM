/**
 * @file BasicNames_test.cpp
 * @short
 * @date 18.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/BasicNames.h"

#include "boost/test/unit_test.hpp"

using std::cout;
using std::endl;

namespace natrium {

BOOST_AUTO_TEST_SUITE(BasicNames_test)


BOOST_AUTO_TEST_CASE(BasicNames_pout_Test){
	cout << "BasicNames_pout_Test..." << endl;

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	BOOST_CHECK( ((bool) mpi_rank ) != perr.is_active());
	BOOST_CHECK( ((bool) mpi_rank ) != pout.is_active());

	cout << "done." << endl;
} /* BasicNames_pout_Test */


BOOST_AUTO_TEST_SUITE_END()


} /* namespace natrium */
