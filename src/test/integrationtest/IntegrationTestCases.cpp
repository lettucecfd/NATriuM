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

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/stencils/Stencil.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q27.h"

#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"
#include "natrium/benchmarks/CouetteFlow2D.h"
#include "natrium/benchmarks/CouetteFlow3D.h"
#include "natrium/benchmarks/AdvectionBenchmark.h"

namespace natrium {
namespace IntegrationTestCases {

TestResult ConvergencePureLinearAdvectionSmooth() {

	TestResult result;
	result.id = 1;
	result.name = "Convergence Test: Pure linear advection (smooth problem)";
	result.details =
			"This test runs the advection solver for a smooth periodic sine profile."
					"Theoretically, the solver has to converge exponentially, which is tested.";
	result.time = clock();
	// smooth: 100*(0.3*dx)**(p+1)
	// nonsmooth: (0.5*p)**(-1.5)*4*dx**2

	bool is_smooth = true;
	for (size_t N = 2; N <= 3; N++) {
		for (size_t orderOfFiniteElement = 4; orderOfFiniteElement <= 8;
				orderOfFiniteElement += 4) {

			double deltaX = 1. / (pow(2, N));
			double deltaT = 0.4 * pow(0.5, N)
					/ ((orderOfFiniteElement + 1) * (orderOfFiniteElement + 1));
			double t_end = 0.1;
			if (t_end/deltaT <= 5) {
				continue;
			}
			AdvectionBenchmark::AdvectionResult advectionResult =
					AdvectionBenchmark::oneTest(N, orderOfFiniteElement, deltaT,
							t_end, RUNGE_KUTTA_5STAGE, NONE, is_smooth, false, false);

			// Analysis
			// Velocity error (compare Paper by Min and Lee)
			std::stringstream stream1;
			stream1 << "|f-f_ref|_sup; N=" << N << "; p=" << orderOfFiniteElement;
			result.quantity.push_back(stream1.str());
			double expected = std::log10(400*std::pow(0.3*deltaX, orderOfFiniteElement+1));
			result.expected.push_back(expected);
			result.threshold.push_back(0.8);
			result.outcome.push_back(std::log10(advectionResult.normSup));

		} /* for p */
	} /* for N */

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


} /* ConvergencePureLinearAdvectionSmooth */

TestResult ConvergencePureLinearAdvectionNonsmooth() {

	TestResult result;
	result.id = 2;
	result.name =
			"Convergence Test: Pure linear advection (non-smooth problem)";
	result.details =
			"This test runs the advection solver for a non-smooth periodic profile, which is "
					"only twice continuously differentiable."
					"The theoretical convergence rate is O(h^min( p+1 , 2 ) * p^(-1.5)) , which is tested.";
	result.time = clock();
	// smooth: 100*(0.3*dx)**(p+1)
	// nonsmooth: (0.5*p)**(-1.5)*4*dx**2

	bool is_smooth = false;
	for (size_t N = 2; N <= 3; N++) {
		for (size_t orderOfFiniteElement = 4; orderOfFiniteElement <= 8;
				orderOfFiniteElement += 4) {

			double deltaX = 1. / (pow(2, N));
			double deltaT = 0.4 * pow(0.5, N)
					/ ((orderOfFiniteElement + 1) * (orderOfFiniteElement + 1));
			double t_end = 0.1;
			if (t_end/deltaT <= 5) {
				continue;
			}
			AdvectionBenchmark::AdvectionResult advectionResult =
					AdvectionBenchmark::oneTest(N, orderOfFiniteElement, deltaT,
							t_end, RUNGE_KUTTA_5STAGE, NONE, is_smooth, false, false);

			// Analysis
			std::stringstream stream1;
			stream1 << "log10(|f-f_ref|_sup); N=" << N << "; p=" << orderOfFiniteElement;
			result.quantity.push_back(stream1.str());
			double expected = std::log10(0.5 * std::pow(0.5*orderOfFiniteElement, -1.5) * deltaX*deltaX);
			result.expected.push_back(expected);
			result.threshold.push_back(0.8);
			result.outcome.push_back( std::log10(advectionResult.normSup));

		} /* for p */
	} /* for N */

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

}
/* ConvergencePureLinearAdvectionNonsmooth */


TestResult ConvergenceTestPeriodic() {

	TestResult result;
	result.id = 3;
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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1. / Ma);
	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getMesh(), orderOfFiniteElement, D2Q9(scaling),
			CFL);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(1.0 / (2 * viscosity) / dt);
	configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
	//configuration->setCollisionScheme(KBC_STANDARD);

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
	result.id = 4;
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
	const double CFL = 5.0;

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1. / Ma);
	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getMesh(), orderOfFiniteElement, D2Q9(scaling),
			CFL);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(1.0 / (2 * viscosity) / dt);
	//configuration->setTimeIntegrator(THETA_METHOD);
	//configuration->setThetaMethodTheta(0.5);
	configuration->setDealIntegrator(SDIRK_TWO_STAGES);
	configuration->setTimeIntegrator(OTHER);
	//configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
	//configuration->setCollisionScheme(KBC_STANDARD);

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

TestResult ConvergenceTestExponentialLBM() {

	TestResult result;
	result.id = 5;
	result.name = "Convergence Test: Exponential time stepping";
	result.details =
			"This test runs the Taylor Green vortex benchmark on a 8x8 grid with FE order 4 and CFL=5.0"
					" It compares the simulated decay of kinetic energy with the analytic solution."
					" The kinematic viscosity is nu=1 and the Reynolds number 2*PI."
					" The simulated and reference E_kin are compared at t=1/(2 nu)";
	result.time = clock();

	// Initialization
	const double viscosity = 1;
	const double Ma = 0.05;
	const double orderOfFiniteElement = 4;
	const double scaling = sqrt(3) * 1 / Ma;
	const double refinementLevel = 3;
	const double CFL = 5.0;

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1. / Ma);
	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getMesh(), orderOfFiniteElement, D2Q9(scaling),
			CFL);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(1.0 / (2 * viscosity) / dt);
	configuration->setTimeIntegrator(EXPONENTIAL);
	//configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
	//configuration->setCollisionScheme(KBC_STANDARD);

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
} /* ConvergenceTestExponentialLBM */

TestResult ConvergenceTestDealIIWrapper() {

	TestResult result;
	result.id = 6;
	result.name = "Convergence Test: DealIIWrapper time stepping (concretely: SDIRK)";
	result.details =
			"This test runs the Taylor Green vortex benchmark on a 8x8 grid with FE order 4 and CFL=5.0"
					" It compares the simulated decay of kinetic energy with the analytic solution."
					" The kinematic viscosity is nu=1 and the Reynolds number 2*PI."
					" The simulated and reference E_kin are compared at t=1/(2 nu)";
	result.time = clock();

	// Initialization
	const double viscosity = 1;
	const double Ma = 0.05;
	const double orderOfFiniteElement = 4;
	const double scaling = sqrt(3) * 1 / Ma;
	const double refinementLevel = 3;
	const double CFL = 5.0;

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1. / Ma);
	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getMesh(), orderOfFiniteElement, D2Q9(scaling),
			CFL);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(1.0 / (2 * viscosity) / dt);
	configuration->setTimeIntegrator(OTHER);
	configuration->setDealIntegrator(SDIRK_TWO_STAGES);
	//configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
	//configuration->setCollisionScheme(KBC_STANDARD);

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
} /* ConvergenceTestDealIIWrapper */

TestResult ConvergenceTest3D() {

	TestResult result;
	result.id = 7;
	result.name = "Convergence Test: 3D flow";
	result.details =
			"This test runs the Couette flow benchmark on a 2x2x2 grid with Re = 10 and  CFL=0.4."
					"The exponential convergence of maximum velocity errors is tested.";
	result.time = clock();

	// Initialization (with standard LB units <=> scaling=1)
	const double Re = 10;
	const double Ma = 0.1;
	double scaling = 1;
	const double U = scaling * Ma / sqrt(3);
	const double L = 1.0;
	const double viscosity = U * L / Re;
	const double refinementLevel = 1;
	const double CFL = 1.0;
	const double t0 = 1.0;

	boost::shared_ptr<Benchmark<3> > benchmark = boost::make_shared<CouetteFlow3D>(viscosity,
			U, refinementLevel, L, t0);

	size_t orderOfFiniteElement = 4;
	// Initialization
	double dt = CFDSolverUtilities::calculateTimestep<3>(
			*benchmark->getMesh(), orderOfFiniteElement,
			D3Q27(scaling), CFL);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setStencil(Stencil_D3Q15);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setTimeStepSize(dt);
	configuration->setTimeIntegrator(RUNGE_KUTTA_5STAGE);
	configuration->setConvergenceThreshold(1e-3);

	// Simulation (simulate 1 time unit from t=40.0)
	BenchmarkCFDSolver<3> solver(configuration, benchmark);
	solver.getSolverStats()->update();
	solver.run();
	solver.getErrorStats()->update();

	// Analysis
	// Velocity: 1 percent allowed
	result.quantity.push_back("|u-u_ref|_2/|u_ref|_2");
	result.expected.push_back(0);
	result.threshold.push_back(1e-2);
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
} /* ConvergenceTest3D */


TestResult ConvergenceTestMovingWall() {
	TestResult result;
	result.id = 8;
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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<CouetteFlow2D>(viscosity,
			U, refinementLevel, L, t0);

	/*for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 12;
	 orderOfFiniteElement += 2) {*/
	for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 10;
			orderOfFiniteElement += 2) {
		// Initialization
		double dt = CFDSolverUtilities::calculateTimestep<2>(
				*benchmark->getMesh(), orderOfFiniteElement,
				D2Q9(scaling), CFL);
		boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
				SolverConfiguration>();
		configuration->setSwitchOutputOff(true);
		//configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setTimeStepSize(dt);
		configuration->setTimeIntegrator(RUNGE_KUTTA_5STAGE);
		configuration->setNumberOfTimeSteps(1.0 / dt);
		//configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
		//configuration->setCollisionScheme(KBC_STANDARD);

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
		result.expected.push_back(
				0.002 * pow(2.0, -(orderOfFiniteElement + 1.0)));
		result.threshold.push_back(
				0.002 * pow(2.0, -(orderOfFiniteElement + 1.0)));
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

