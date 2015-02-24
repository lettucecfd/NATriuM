/*
 * IntegrationTestCases.cpp
 *
 *  Created on: Jan 28, 2015
 *      Author: kraemer
 */

#include "IntegrationTestCases.h"

#include <time.h>
#include <math.h>
#include <sstream>

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/CFDSolverUtilities.h"

#include "step-1/TaylorGreenVortex2D.h"
#include "step-2/CouetteFlow2D.h"

namespace natrium {
namespace IntegrationTestCases {

TestResult ConvergenceTestPeriodic() {

	TestResult result;
	result.id = 1;
	result.name = "Convergence Test: Periodic Boundaries";
	result.details =
			"This test runs the Taylor Green vortex benchmark on a 8x8 grid with FE order 4 and CFL=0.4."
					"It compares the simulated decay of kinetic energy with the analytic solution."
					"The kinematic viscosity is nu=1 and the Reynolds number 2*PI."
					"The simulated and reference E_kin are compared at t=1/(2 nu)";
	result.time = clock();

	// Initialization
	const double viscosity = 1;
	const double Ma = 0.05;
	const double orderOfFiniteElement = 4;
	const double scaling = sqrt(3) * 1 / Ma;
	const double refinementLevel = 3;
	const double CFL = 0.4;

	shared_ptr<Benchmark<2> > benchmark = make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1./Ma);
	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getTriangulation(), orderOfFiniteElement,
			D2Q9IncompressibleModel(scaling), CFL);

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(1.0 / (2 * viscosity) / dt);

	// Simulation
	BenchmarkCFDSolver<2> solver(configuration, benchmark);
	solver.getSolverStats()->update();
	double Ekin0 = solver.getSolverStats()->getKinE();
	solver.run();
	solver.getSolverStats()->update();
	solver.getErrorStats()->update();

	// Analysis
	// Kinetic Energy
	result.quantity.push_back("E_kin(t=1/(2 nu))/E_kin(t=0)");
	result.expected.push_back(exp(-2));
	result.threshold.push_back(1e-2);
	result.outcome.push_back(solver.getSolverStats()->getKinE() / Ekin0);

	// Velocity: 5 percent allowed
	result.quantity.push_back("|u-u_ref|_2/|u_ref|_2");
	result.expected.push_back(0);
	result.threshold.push_back(5e-2);
	result.outcome.push_back(
			solver.getErrorStats()->getL2VelocityError()
					/ solver.getErrorStats()->getL2UAnalytic());

	// Finalize test
	result.time = (clock() - result.time) / CLOCKS_PER_SEC;
	assert(result.quantity.size() == result.expected.size());
	assert(result.quantity.size() == result.threshold.size());
	assert(result.quantity.size() == result.outcome.size());
	result.success = true;
	for (size_t i = 0; i < result.quantity.size(); i++) {
		if (fabs(result.expected.at(i) - result.outcome.at(i))
				> result.threshold.at(i)) {
			result.success = false;
			*result.error_msg << result.quantity.at(i)
					<< " not below threshold.";
		}
	}
	return result;
} /* ConvergenceTestPeriodicBoundary */


TestResult ConvergenceTestImplicitLBM() {

	TestResult result;
	result.id = 2;
	result.name = "Convergence Test: Implicit time stepping";
	result.details =
			"This test runs the Taylor Green vortex benchmark on a 8x8 grid with FE order 4 and CFL=5."
					" It compares the simulated decay of kinetic energy with the analytic solution."
					" The particular time stepping scheme is a theta method with theta=0.5."
					" The kinematic viscosity is nu=1 and the Reynolds number 2*PI."
					" The simulated and reference E_kin are compared at t=1/(2 nu)";
	result.time = clock();

	// Initialization
	const double viscosity = 1;
	const double Ma = 0.05;
	const double orderOfFiniteElement = 4;
	const double scaling = sqrt(3) * 1 / Ma;
	const double refinementLevel = 3;
	const double CFL = 1.0;

	shared_ptr<Benchmark<2> > benchmark = make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1./Ma);
	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getTriangulation(), orderOfFiniteElement,
			D2Q9IncompressibleModel(scaling), CFL);

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(1.0 / (2 * viscosity) / dt);
	//configuration->setTimeIntegrator(THETA_METHOD);
	//configuration->setThetaMethodTheta(0.5);
	configuration->setDealIntegrator(SDIRK_TWO_STAGES);
	configuration->setTimeIntegrator(OTHER);

	// Simulation
	BenchmarkCFDSolver<2> solver(configuration, benchmark);
	solver.getSolverStats()->update();
	double Ekin0 = solver.getSolverStats()->getKinE();
	solver.run();
	solver.getSolverStats()->update();
	solver.getErrorStats()->update();

	// Analysis
	// Kinetic Energy
	result.quantity.push_back("E_kin(t=1/(2 nu))/E_kin(t=0)");
	result.expected.push_back(exp(-2));
	result.threshold.push_back(1e-2);
	result.outcome.push_back(solver.getSolverStats()->getKinE() / Ekin0);

	// Velocity: 5 percent allowed
	result.quantity.push_back("|u-u_ref|_2/|u_ref|_2");
	result.expected.push_back(0);
	result.threshold.push_back(5e-2);
	result.outcome.push_back(
			solver.getErrorStats()->getL2VelocityError()
					/ solver.getErrorStats()->getL2UAnalytic());

	// Finalize test
	result.time = (clock() - result.time) / CLOCKS_PER_SEC;
	assert(result.quantity.size() == result.expected.size());
	assert(result.quantity.size() == result.threshold.size());
	assert(result.quantity.size() == result.outcome.size());
	result.success = true;
	for (size_t i = 0; i < result.quantity.size(); i++) {
		if (fabs(result.expected.at(i) - result.outcome.at(i))
				> result.threshold.at(i)) {
			result.success = false;
			*result.error_msg << result.quantity.at(i)
					<< " not below threshold.";
		}
	}
	return result;
} /* ConvergenceTestImplicitLBM */


TestResult ConvergenceTestMovingWall() {
	TestResult result;
	result.id = 3;
	result.name = "Convergence Test: Wall Boundaries";
	result.details =
			"This test runs the Unsteady Couette flow benchmark on a 4x4 grid with FE orders 2,4,6,8,10,12 "
					"and CFL=0.4. It reproduces exponential convergence observed by Min and Lee for Re=2000."
					"The stencil scaling, Mach and Reynolds numbers are 1, 0.05 and 2000, respectively."
					"The simulation STARTS AT t=40.0 and lasts only one time unit (!!!)"
	 	 	 	 	"The exponential convergence of maximum velocity errors is tested.";
	result.time = clock();

	// Initialization (with standard LB units <=> scaling=1)
	const double Re = 2000;
	const double Ma = 0.05;
	double scaling = 1;
	const double U = scaling * Ma / sqrt(3);
	const double L = 1.0;
	const double viscosity = U * L / Re;
	const double refinementLevel = 2;
	const double CFL = 0.4;
	const double t0 = 40.0;

	shared_ptr<Benchmark<2> > benchmark = make_shared<CouetteFlow2D>(viscosity,
			U, refinementLevel, L, t0);

	/*for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 12;
			orderOfFiniteElement += 2) {*/
	for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 10;
			orderOfFiniteElement += 2) {
		// Initialization
		double dt = CFDSolverUtilities::calculateTimestep<2>(
				*benchmark->getTriangulation(), orderOfFiniteElement,
				D2Q9IncompressibleModel(scaling), CFL);
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		configuration->setSwitchOutputOff(true);
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setTimeStepSize(dt);
		configuration->setNumberOfTimeSteps(1.0 / dt);

		// Simulation (simulate 1 time unit from t=40.0)
		BenchmarkCFDSolver<2> solver(configuration, benchmark);
		solver.getSolverStats()->update();
		solver.run();
		solver.getErrorStats()->update();

		// Analysis
		// Velocity error (compare Paper by Min and Lee)
		std::stringstream stream1;
		stream1 << "|u-u_ref|_sup; p=" << orderOfFiniteElement;
		result.quantity.push_back(stream1.str());
		result.expected.push_back(0.002 * pow(2.0, -( orderOfFiniteElement + 1.0)));
		result.threshold.push_back(0.002 * pow(2.0, -( orderOfFiniteElement + 1.0)));
		result.outcome.push_back(solver.getErrorStats()->getMaxVelocityError());
	}

	// Finalize test
	result.time = (clock() - result.time) / CLOCKS_PER_SEC;
	assert(result.quantity.size() == result.expected.size());
	assert(result.quantity.size() == result.threshold.size());
	assert(result.quantity.size() == result.outcome.size());
	result.success = true;
	for (size_t i = 0; i < result.quantity.size(); i++) {
		if (fabs(result.expected.at(i) - result.outcome.at(i))
				> result.threshold.at(i)) {
			result.success = false;
			*result.error_msg << result.quantity.at(i)
					<< " not below threshold.";
		}
	}
	return result;
} /* ConvergenceTestMovingWall */


} /* namespace IntegrationTests */
} /* namespace natrium */

