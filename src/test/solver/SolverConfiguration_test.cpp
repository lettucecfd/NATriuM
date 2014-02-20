/**
 * @file SolverConfiguration_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "solver/SolverConfiguration.h"

#include "boost/test/unit_test.hpp"

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(CFDSolverConfiguration_test)

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_Construction_test){
	cout << "CFDSolverConfiguration_Construction_test..." << endl;
	SolverConfiguration solver;
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_OutputFlgas_test){
	cout << "CFDSolverConfiguration_OutputFlgas_test..." << endl;
	SolverConfiguration cfg;

	cout << "done";
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
