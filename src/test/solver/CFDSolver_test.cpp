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
	SteadyPeriodicTestFlow2D testFlow(tau);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_Construction_test) {
	cout << "CFDSolver_Construction_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	double tau = 0.9;
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<SteadyPeriodicTestFlow2D>(tau);
	CFDSolver<2> solver(testConfiguration, testFlow);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_SteadyStreaming_test) {
	cout << "CFDSolver_SteadyStreaming_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	double tau = 0.9;
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<SteadyPeriodicTestFlow2D>(tau);
	CFDSolver<2> solver(testConfiguration, testFlow);
	solver.run();
	// check results (must be the same as before)
	const distributed_vector& rho = solver.getDensity();
	const vector<distributed_vector>& v = solver.getVelocity();
	for (size_t i = 0; i < rho.size(); i++){
		BOOST_CHECK(fabs(rho(i) - 1.0) < 1e-5);
		BOOST_CHECK(fabs(v.at(0)(i) - 0.1) < 1e-5);
		BOOST_CHECK(fabs(v.at(1)(i) - 0.1) < 1e-5);
	}
	cout << "done" << endl;
}


BOOST_AUTO_TEST_CASE(CFDSolver_UnsteadyStreaming_test) {
	cout << "CFDSolver_UnsteadyStreaming_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	double tau = 0.9;

	shared_ptr<ProblemDescription<2> > testFlow = make_shared<UnsteadyPeriodicTestFlow2D>(tau);
	CFDSolver<2> solver(testConfiguration, testFlow);
	solver.run();
	// check results (must be the same as before)
	const distributed_vector& rho = solver.getDensity();
	const vector<distributed_vector>& v = solver.getVelocity();
	for (size_t i = 0; i < rho.size(); i++){
		//BOOST_CHECK(fabs(rho(i) - 1.0) < 1e-5);
		//BOOST_CHECK(fabs(v.at(0)(i) - 0.1) < 1e-5);
		//BOOST_CHECK(fabs(v.at(1)(i) - 0.1) < 1e-5);
	}
	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
