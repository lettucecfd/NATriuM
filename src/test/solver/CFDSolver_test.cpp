/**
 * @file CFDSolver_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/solver/CFDSolver.h"

#include "boost/test/included/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"

#include "PeriodicTestFlow2D.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(CFDSolver_test)

BOOST_AUTO_TEST_CASE(CFDSolver_CreateTestFlow_test) {
	pout << "CFDSolver_CreateTestFlow_test..." << endl;
	double viscosity = 0.9;
	size_t refinementLevel = 3;
	SteadyPeriodicTestFlow2D testFlow(viscosity, refinementLevel);
	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_Construction_test) {
	pout << "CFDSolver_Construction_test..." << endl;
	boost::shared_ptr<SolverConfiguration> testConfiguration =
			boost::make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	testConfiguration->setCommandLineVerbosity(DETAILED);
	testConfiguration->setUserInteraction(false);
	size_t refinementLevel = 3;
	double viscosity = 0.9;
	boost::shared_ptr<ProblemDescription<2> > testFlow = boost::make_shared<
			SteadyPeriodicTestFlow2D>(viscosity, refinementLevel);
	CFDSolver<2> solver(testConfiguration, testFlow);
	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_SteadyStreaming_test) {
	pout << "CFDSolver_SteadyStreaming_test..." << endl;
	boost::shared_ptr<SolverConfiguration> testConfiguration =
			boost::make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	size_t refinementLevel = 3;
	double deltaX = 1.
			/ (pow(2, refinementLevel)
					* (testConfiguration->getSedgOrderOfFiniteElement()));
	testConfiguration->setCFL(0.4);
	testConfiguration->setNumberOfTimeSteps(100);
	// set viscosity so that tau = 1
	double viscosity = 0.5 * deltaX / 3;
	boost::shared_ptr<ProblemDescription<2> > testFlow = boost::make_shared<
			SteadyPeriodicTestFlow2D>(viscosity, refinementLevel);

	CFDSolver<2> solver(testConfiguration, testFlow);
	solver.run();
	// check results (must be the same as before)
	const distributed_vector& rho = solver.getDensity();
	const vector<distributed_vector>& v = solver.getVelocity();

	//for all degrees of freedom on current processor
	const dealii::IndexSet& locally_owned_dofs =
			solver.getAdvectionOperator()->getLocallyOwnedDofs();
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;
		BOOST_CHECK(fabs(rho(i) - 1.0) < 1e-5);
		BOOST_CHECK(fabs(v.at(0)(i) - 0.1) < 1e-5);
		BOOST_CHECK(fabs(v.at(1)(i) - 0.1) < 1e-5);
	}
	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_UnsteadyStreaming_test) {
	pout << "CFDSolver_UnsteadyStreaming_test..." << endl;
	boost::shared_ptr<SolverConfiguration> testConfiguration =
			boost::make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	size_t refinementLevel = 3;
	testConfiguration->setCFL(0.8);
	testConfiguration->setNumberOfTimeSteps(200);
	// set viscosity so that tau = 1
	double viscosity = 1. / 5;
	testConfiguration->setStencilScaling(1);

	boost::shared_ptr<ProblemDescription<2> > testFlow = boost::make_shared<
			UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);

	CFDSolver<2> solver(testConfiguration, testFlow);
	solver.run();
	// check results (must be nearly at equilibrium)
	const distributed_vector& rho = solver.getDensity();
	const vector<distributed_vector>& v = solver.getVelocity();
	//for all degrees of freedom on current processor
	const dealii::IndexSet& locally_owned_dofs =
			solver.getAdvectionOperator()->getLocallyOwnedDofs();
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;
		BOOST_CHECK(fabs(rho(i) - 1.0) < 1e-5);
		BOOST_CHECK(fabs(v.at(0)(i) - 0.0) < 0.001);
		BOOST_CHECK(fabs(v.at(1)(i) - 0.0) < 0.001);
	}
	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolver_Restart_test) {
	pout << "CFDSolver_Restart_test..." << endl;

	// Solver configuration
	string directory = "/tmp/test-restart";
	boost::shared_ptr<SolverConfiguration> testConfiguration =
			boost::make_shared<SolverConfiguration>();
	testConfiguration->setOutputDirectory(directory);
	testConfiguration->setOutputCheckpointInterval(10);
	size_t refinementLevel = 3;
	testConfiguration->setCFL(0.1);
	testConfiguration->setNumberOfTimeSteps(15);
	double viscosity = 1. / 5;
	testConfiguration->setStencilScaling(1);

	// create problem and solver solver
	boost::shared_ptr<ProblemDescription<2> > testFlow = boost::make_shared<
			UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setCommandLineVerbosity(0);
	CFDSolver<2> solver(testConfiguration, testFlow);

	// first run
	solver.run();

	// restart run
	boost::shared_ptr<ProblemDescription<2> > testFlow2 = boost::make_shared<
				UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);
	testConfiguration->setRestartAtIteration(10);
	testConfiguration->setNumberOfTimeSteps(25);
	CFDSolver<2> solver2(testConfiguration, testFlow2);
	BOOST_CHECK(solver2.getIterationStart() == 10);
	solver2.run();

	pout << "done" << endl;
} // CFDSolver_Restart_test

BOOST_AUTO_TEST_CASE(CFDSolver_IterativeInit_test) {
	pout << "CFDSolver_IterativeInit_test..." << endl;
	boost::shared_ptr<SolverConfiguration> testConfiguration =
			boost::make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	testConfiguration->setInitializationScheme(ITERATIVE);
	testConfiguration->setIterativeInitializationNumberOfIterations(10);
	testConfiguration->setIterativeInitializationResidual(0.001);
	size_t refinementLevel = 3;
	testConfiguration->setCFL(0.4);
	testConfiguration->setNumberOfTimeSteps(100);
	// set viscosity so that tau = 1
	double viscosity = 1. / 5;
	testConfiguration->setStencilScaling(1);

	boost::shared_ptr<ProblemDescription<2> > testFlow = boost::make_shared<
			UnsteadyPeriodicTestFlow2D>(viscosity, refinementLevel);

	CFDSolver<2> solver(testConfiguration, testFlow);

	pout << "done" << endl;
} // CFDSolver_IterativeInit_test

BOOST_AUTO_TEST_CASE(CFDSolver_StopCondition_test) {
	pout << "CFDSolver_StopCondition_test..." << endl;

	boost::shared_ptr<SolverConfiguration> testConfiguration =
			boost::make_shared<SolverConfiguration>();
	testConfiguration->setSwitchOutputOff(true);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setSimulationEndTime(5);
	testConfiguration->setNumberOfTimeSteps(10);
	testConfiguration->setConvergenceThreshold(5e-3);

	size_t refinementLevel = 2;
	testConfiguration->setCFL(0.1);
	// set viscosity so that tau = 1
	double viscosity = 1. / 5;

	boost::shared_ptr<ProblemDescription<2> > testFlow = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, refinementLevel);
	CFDSolver<2> solver(testConfiguration, testFlow);

	solver.run();
	BOOST_CHECK_EQUAL(solver.getIteration(), size_t(10));
	testConfiguration->setNumberOfTimeSteps(1000000);

	solver.run();
	BOOST_CHECK_GE(solver.getTime(), 5);
	BOOST_CHECK_LT(solver.getTime(), 5 + solver.getTimeStepSize());
	testConfiguration->setSimulationEndTime(1000000.0);

	solver.run();
	//BOOST_CHECK_LE(solver.getResiduumDensity(), 5e-2);
	BOOST_CHECK_LE(solver.getResiduumVelocity(), 5e-2);

	pout << "done" << endl;
} /* CFDSolver_StopCondition_test */

BOOST_AUTO_TEST_SUITE_END()

