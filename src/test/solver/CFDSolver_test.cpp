/**
 * @file CFDSolver_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "solver/CFDSolver.h"

#include "boost/test/unit_test.hpp"

#include "utilities/BasicNames.h"

#include "PeriodicTestFlow2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(CFDSolver_test)

BOOST_AUTO_TEST_CASE(CFDSolver_CreateTestFlow_test) {
	cout << "CFDSolver_CreateTestFlow_test..." << endl;
	double tau = 0.9;
	numeric_vector u(2);
	u(0) = 0.5;
	u(1) = 0.5;
	PeriodicTestFlow2D testFlow(tau, u);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_Construction_test) {
	cout << "CFDSolver_Construction_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	double tau = 0.9;
	numeric_vector u(2);
	u(0) = 0.5;
	u(1) = 0.5;
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<PeriodicTestFlow2D>(tau, u);
	CFDSolver<2> solver(testConfiguration, testFlow);
	cout << "done" << endl;
}



BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
