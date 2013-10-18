/**
 * @file BasicNames_test.cpp
 * @short
 * @date 18.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "utilities/BasicNames.h"

#include "boost/test/unit_test.hpp"

using std::cout;
using std::endl;

namespace natrium {

BOOST_AUTO_TEST_SUITE(BasicNames_test)


#ifdef WITH_PETSC
BOOST_AUTO_TEST_CASE(BasicNames_ParallelPETSc_Test){
	cout << "BasicNames_ParallelPETSc_Test..." << endl;

	distributed_vector vector;

	cout << "done..." << endl;
} /* BasicNames_ParallelPETSc_Test */
#endif


BOOST_AUTO_TEST_SUITE_END()


} /* namespace natrium */
