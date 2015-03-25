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
#include "TaylorGreenVortex2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(CFDSolver_test)

BOOST_AUTO_TEST_CASE(CFDSolver_CreateTestFlow_test) {
	cout << "CFDSolver_CreateTestFlow_test..." << endl;
	double viscosity = 0.9;
	size_t refinementLevel = 3;
	SteadyPeriodicTestFlow2D testFlow(viscosity, refinementLevel);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_Construction_test) {
	cout << "CFDSolver_Construction_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	size_t refinementLevel = 3;
	double viscosity = 0.9;
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<SteadyPeriodicTestFlow2D>(viscosity, refinementLevel);
	CFDSolver<2> solver(testConfiguration, testFlow);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_SteadyStreaming_test) {
	cout << "CFDSolver_SteadyStreaming_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	size_t refinementLevel = 3;
	double deltaX = 1./(pow(2,refinementLevel)*(testConfiguration->getSedgOrderOfFiniteElement()));
	testConfiguration->setTimeStepSize(0.5*deltaX);
	testConfiguration->setNumberOfTimeSteps(100);
	// set viscosity so that tau = 1
	double viscosity = 0.5*deltaX/3;
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<SteadyPeriodicTestFlow2D>(viscosity, refinementLevel);

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
	testConfiguration->setSwitchOutputOff(true);
	size_t refinementLevel = 3;
	double deltaX = 1./(pow(2,refinementLevel)*(testConfiguration->getSedgOrderOfFiniteElement()));
	testConfiguration->setTimeStepSize(0.1*deltaX);
	testConfiguration->setNumberOfTimeSteps(100);
	// set viscosity so that tau = 1
	double viscosity = 1./5;
	testConfiguration->setStencilScaling(sqrt(3*viscosity/(testConfiguration->getTimeStepSize())));

	shared_ptr<ProblemDescription<2> > testFlow = make_shared<UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);

	CFDSolver<2> solver(testConfiguration, testFlow);
	solver.run();
	// check results (must be nearly at equilibrium)
	const distributed_vector& rho = solver.getDensity();
	const vector<distributed_vector>& v = solver.getVelocity();
	for (size_t i = 0; i < rho.size(); i++){
		BOOST_CHECK(fabs(rho(i) - 1.0) < 1e-5);
		BOOST_CHECK(fabs(v.at(0)(i) - 0.0) < 0.001);
		BOOST_CHECK(fabs(v.at(1)(i) - 0.0) < 0.001);
	}
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_Restart_test) {
	cout << "CFDSolver_Restart_test..." << endl;

	// Solver configuration
	string directory = "/tmp/test-restart";
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	testConfiguration->setOutputDirectory(directory);
	testConfiguration->setOutputCheckpointInterval(10);
	size_t refinementLevel = 3;
	double deltaX = 1./(pow(2,refinementLevel)*(testConfiguration->getSedgOrderOfFiniteElement()));
	testConfiguration->setTimeStepSize(0.1*deltaX);
	testConfiguration->setNumberOfTimeSteps(15);
	double viscosity = 1./5;
	testConfiguration->setStencilScaling(sqrt(3*viscosity/(testConfiguration->getTimeStepSize())));

	// create problem and solver solver
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setCommandLineVerbosity(0);
	CFDSolver<2> solver(testConfiguration, testFlow);

	// first run
	solver.run();

	// restart run
	testConfiguration->setRestartAtLastCheckpoint(true);
	testConfiguration->setNumberOfTimeSteps(25);
	CFDSolver<2> solver2(testConfiguration, testFlow);
	BOOST_CHECK(solver2.getIterationStart() == 10);
	solver2.run();


	cout << "done" << endl;
} /* CFDSolver_Restart_test */

BOOST_AUTO_TEST_CASE(CFDSolver_IterativeInit_test) {
	cout << "CFDSolver_IterativeInit_test..." << endl;
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	testConfiguration->setInitializationScheme(ITERATIVE);
	testConfiguration->setIterativeInitializationNumberOfIterations(10);
	testConfiguration->setIterativeInitializationResidual(0.001);
	size_t refinementLevel = 3;
	double deltaX = 1./(pow(2,refinementLevel)*(testConfiguration->getSedgOrderOfFiniteElement()));
	testConfiguration->setTimeStepSize(0.1*deltaX);
	testConfiguration->setNumberOfTimeSteps(100);
	// set viscosity so that tau = 1
	double viscosity = 1./5;
	testConfiguration->setStencilScaling(sqrt(3*viscosity/(testConfiguration->getTimeStepSize())));

	shared_ptr<ProblemDescription<2> > testFlow = make_shared<UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);

	CFDSolver<2> solver(testConfiguration, testFlow);

	cout << "done" << endl;
} /* CFDSolver_IterativeInit_test */

BOOST_AUTO_TEST_CASE(CFDSolver_StopCondition_test) {
	cout << "CFDSolver_StopCondition_test..." << endl;

	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(false);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setSimulationEndTime(0.5);
	testConfiguration->setNumberOfTimeSteps(10);
	testConfiguration->setConvergenceThreshold(5e-2);

	size_t refinementLevel = 2;
	double deltaX = 1./(pow(2,refinementLevel)*(testConfiguration->getSedgOrderOfFiniteElement()));
	testConfiguration->setTimeStepSize(0.1*deltaX);
	// set viscosity so that tau = 1
	double viscosity = 1./5;

	shared_ptr<ProblemDescription<2> > testFlow = make_shared<TaylorGreenVortex2D>(viscosity, refinementLevel);
	CFDSolver<2> solver(testConfiguration, testFlow);

	solver.run();
	BOOST_CHECK(solver.getIteration() == 10);
	testConfiguration->setNumberOfTimeSteps(1000000);

	solver.run();
	BOOST_CHECK(solver.getTime() >= 0.5);
	BOOST_CHECK(solver.getTime() < 0.5 + 0.1*deltaX);
	testConfiguration->setSimulationEndTime(1000000.0);

	solver.run();
	BOOST_CHECK_LE(solver.getResiduumDensity(), 5e-2);
	BOOST_CHECK_LE(solver.getResiduumVelocity(), 5e-2);


	cout << "done" << endl;
} /* CFDSolver_StopCondition_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
