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
#include "natrium/benchmarks/PoiseuilleFlow2D.h"
#include "natrium/benchmarks/PoiseuilleFlow3D.h"

#include "natrium/benchmarks/AdvectionBenchmark.h"


using namespace natrium;
namespace IntegrationTestCases {

TestResult ConvergenceSEDGLinearAdvectionSmooth() {

	TestResult result;
	result.id = 1;
	result.name = "Convergence Test: SEDG linear advection (smooth problem)";
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
			if (t_end / deltaT <= 5) {
				continue;
			}
			AdvectionBenchmark::AdvectionResult advectionResult =
					AdvectionBenchmark::oneTest(N, orderOfFiniteElement, deltaT,
							t_end, RUNGE_KUTTA_5STAGE, NONE, is_smooth, false,
							false);

			// Analysis
			// Velocity error (compare Paper by Min and Lee)
			std::stringstream stream1;
			stream1 << "|f-f_ref|_sup; N=" << N << "; p="
					<< orderOfFiniteElement;
			result.quantity.push_back(stream1.str());
			double expected = std::log10(
					400 * std::pow(0.3 * deltaX, orderOfFiniteElement + 1));
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

} /* ConvergenceSEDGLinearAdvectionSmooth */

TestResult ConvergenceSEDGLinearAdvectionNonsmooth() {

	TestResult result;
	result.id = 2;
	result.name =
			"Convergence Test: SEDG linear advection (non-smooth problem)";
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
			if (t_end / deltaT <= 5) {
				continue;
			}
			AdvectionBenchmark::AdvectionResult advectionResult =
					AdvectionBenchmark::oneTest(N, orderOfFiniteElement, deltaT,
							t_end, RUNGE_KUTTA_5STAGE, NONE, is_smooth, false,
							false);

			// Analysis
			std::stringstream stream1;
			stream1 << "log10(|f-f_ref|_sup); N=" << N << "; p="
					<< orderOfFiniteElement;
			result.quantity.push_back(stream1.str());
			double expected = std::log10(
					0.5 * std::pow(0.5 * orderOfFiniteElement, -1.5) * deltaX
							* deltaX);
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

}
/* ConvergenceSEDGLinearAdvectionNonsmooth */

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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, refinementLevel, 1. / Ma);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(1.0 / (2 * viscosity));
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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, refinementLevel, 1. / Ma);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(1.0 / (2 * viscosity));
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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, refinementLevel, 1. / Ma);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(1.0 / (2 * viscosity));
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
	result.name =
			"Convergence Test: DealIIWrapper time stepping (concretely: SDIRK)";
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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, refinementLevel, 1. / Ma);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(1.0 / (2 * viscosity));
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
	const double refinementLevel = 2;
	const double CFL = 1.0;
	const double t0 = 1.0;

	boost::shared_ptr<Benchmark<3> > benchmark = boost::make_shared<
			CouetteFlow3D>(viscosity, U, refinementLevel, L, t0);

	size_t orderOfFiniteElement = 2;
	// Initialization
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setCommandLineVerbosity(ALL);
	configuration->setStencil(Stencil_D3Q15);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
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

	/*for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 12;
	 orderOfFiniteElement += 2) {*/
	for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 10;
			orderOfFiniteElement += 2) {
		// Initialization
		boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<
				CouetteFlow2D>(viscosity, U, refinementLevel, L, t0);
		boost::shared_ptr<SolverConfiguration> configuration =
				boost::make_shared<SolverConfiguration>();
		configuration->setSwitchOutputOff(true);
		//configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCFL(CFL);
		configuration->setTimeIntegrator(RUNGE_KUTTA_5STAGE);
		configuration->setSimulationEndTime(1.0);
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

TestResult ConvergenceTestForcingSchemes2D() {
	TestResult result;
	result.id = 9;
	result.name = "Convergence Test: Forcing Schemes in 2D";
	result.details =
			"This test runs the Poiseuille flow benchmark on a 2x2 grid with FE order 2."
					"The flow is initialized with a zero profile and driven by an external force."
					"The maximal velocity is then compared with the theoretical value."
					"The stencil is a D2Q9 stencil."
					"Three different Forcing schemes are tested: Shifting velocity, exact difference, and Guo.";
	result.time = clock();

	// Initialization (with standard LB units <=> scaling=1)

	// setup test case
	const double CFL = 0.8;
	const double Re = 10;
	const double u_bulk = 0.0001 / 1.5; //1.0;
	const double height = 3.0;
	const double length = 2.0;
	const double orderOfFiniteElement = 2;
	const double Ma = 0.1;
	const double refinement_level = 1;
	bool is_periodic = true;

	/// create CFD problem
	double viscosity = u_bulk * height / Re;
	const double scaling = sqrt(3) * 1.5 * u_bulk / Ma;

	/// setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setConvergenceThreshold(1e-8);
	//configuration->setTimeIntegrator(OTHER);
	//configuration->setDealIntegrator(SDIRK_TWO_STAGES);

	// SHIFTING VELOCITY
	{
		//pout << "... shifting velocity scheme." << endl;
		configuration->setForcingScheme(SHIFTING_VELOCITY);

		// make solver object and run simulation
		boost::shared_ptr<ProblemDescription<2> > poiseuille2D =
				boost::make_shared<PoiseuilleFlow2D>(viscosity,
						refinement_level, u_bulk, height, length, is_periodic);

		CFDSolver<2> solver(configuration, poiseuille2D);
		solver.run();

		result.quantity.push_back("Shifting Velocity");
		result.expected.push_back(1.5 * u_bulk);
		result.threshold.push_back(0.05 * 1.5 * u_bulk);
		result.outcome.push_back(solver.getMaxVelocityNorm());
		//BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
		//pout << "... done." << endl;

	}

	// EXACT_DIFFERENCE
	{
		//pout << "... exact difference scheme." << endl;
		configuration->setForcingScheme(EXACT_DIFFERENCE);

		// make solver object and run simulation
		boost::shared_ptr<ProblemDescription<2> > poiseuille2D =
				boost::make_shared<PoiseuilleFlow2D>(viscosity,
						refinement_level, u_bulk, height, length, is_periodic);

		CFDSolver<2> solver(configuration, poiseuille2D);
		solver.run();

		result.quantity.push_back("Exact Difference");
		result.expected.push_back(1.5 * u_bulk);
		result.threshold.push_back(0.05 * 1.5 * u_bulk);
		result.outcome.push_back(solver.getMaxVelocityNorm());
		//BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
		//pout << "... done." << endl;
	}

	//GUO
	{
		//pout << "... Guo scheme." << endl;
		configuration->setForcingScheme(GUO);

		// make solver object and run simulation
		boost::shared_ptr<ProblemDescription<2> > poiseuille2D =
				boost::make_shared<PoiseuilleFlow2D>(viscosity,
						refinement_level, u_bulk, height, length, is_periodic);
		CFDSolver<2> solver(configuration, poiseuille2D);
		solver.run();

		result.quantity.push_back("Guo");
		result.expected.push_back(1.5 * u_bulk);
		result.threshold.push_back(0.05 * 1.5 * u_bulk);
		result.outcome.push_back(solver.getMaxVelocityNorm());
		//BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
		//pout << "... done." << endl;
	}

	// Analysis
	// Velocity error (compare Paper by Min and Lee)



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

TestResult ConvergenceTestForcingSchemes3D() {
	TestResult result;
	result.id = 10;
	result.name = "Convergence Test: Forcing Schemes in 3D";
	result.details =
			"This test runs the Poiseuille flow benchmark on a 2x2x2 grid with FE order 2."
					"The flow is initialized with a zero profile and driven by an external force."
					"The maximal velocity is then compared with the theoretical value."
					"Three different Forcing schemes are tested: Shifting velocity, exact difference, and Guo."
					"Two different stencils are used: D3Q19 and D3Q15.";

	result.time = clock();


	// setup test case
	const double CFL = 0.8;
	const double Re = 10;
	const double u_bulk = 0.0001 / 1.5; //1.0;
	const double height = 3.0;
	const double width = 2.0;
	const double length = 2.0;
	const double orderOfFiniteElement = 2;
	const double Ma = 0.1;
	const double refinement_level = 1;
	bool is_periodic = true;

	/// create CFD problem
	double viscosity = u_bulk * height / Re;
	const double scaling = sqrt(3) * 1.5 * u_bulk / Ma;

	/// setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setConvergenceThreshold(1e-8);
	//configuration->setTimeIntegrator(OTHER);
	//configuration->setDealIntegrator(SDIRK_TWO_STAGES);

	// forall stencil types and forcing schemes
	// TODO there is a bug in D3Q27; set "i < 4" when resolved
	for (size_t i = 1; i < 3; i++) {
		configuration->setStencil(static_cast<StencilType>(i));
		for (size_t j = 1; j < 4; j++) {
			configuration->setForcingScheme(static_cast<ForceType>(j));
			std::stringstream s;
			switch (i) {
			case 1: {
				s << "D3Q19; ";
				break;
			}
			case 2: {
				s << "D3Q15; ";
				break;
			}
			case 3: {
				s << "D3Q27; ";
				break;
			}
			}
			switch (j) {
			case 1: {
				s << " shifting velocity scheme; ";
				break;
			}
			case 2: {
				s << " exact difference scheme; ";
				break;
			}
			case 3: {
				s << " Guo scheme; ";
				break;
			}
			}

			// make solver object and run simulation
			boost::shared_ptr<ProblemDescription<3> > poiseuille3D =
					boost::make_shared<PoiseuilleFlow3D>(viscosity,
							refinement_level, u_bulk, height, width, length,
							is_periodic);

			CFDSolver<3> solver(configuration, poiseuille3D);
			solver.run();

			result.quantity.push_back(s.str());
			result.expected.push_back(1.5 * u_bulk);
			result.threshold.push_back(0.05 * 1.5 * u_bulk);
			result.outcome.push_back(solver.getMaxVelocityNorm());

		}
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
}

TestResult ConvergenceTestSemiLagrangianPeriodic() {

	TestResult result;
	result.id = 11;
	result.name = "Convergence Test: Semi-Lagrangian LBM with Periodic Boundaries";
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

	boost::shared_ptr<Benchmark<2> > benchmark = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, refinementLevel, 1. / Ma);

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(1.0 / (2 * viscosity));
	configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
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
} /* ConvergenceTestSemiLagrangianPeriodic */

TestResult ConvergenceTestSemiLagrangianAdvectionSmooth (){

	TestResult result;
	result.id = 12;
	result.name = "Convergence Test: Semi-Lagrangian linear advection (smooth problem)";
	result.details =
			"This test runs the semi-Lagrangian advection solver for a smooth periodic sine profile.";
	result.time = clock();
	// smooth: 100*(0.3*dx)**(p+1)
	// nonsmooth: (0.5*p)**(-1.5)*4*dx**2

	bool is_smooth = true;
	bool semi_lagrangian = true;
	for (size_t N = 2; N <= 3; N++) {
		for (size_t orderOfFiniteElement = 4; orderOfFiniteElement <= 8;
				orderOfFiniteElement += 4) {

			double deltaX = 1. / (pow(2, N));
			double deltaT = 0.4 * pow(0.5, N)
					/ ((orderOfFiniteElement + 1) * (orderOfFiniteElement + 1));
			double t_end = 0.1;
			if (t_end / deltaT <= 5) {
				continue;
			}
			AdvectionBenchmark::AdvectionResult advectionResult =
					AdvectionBenchmark::oneTest(N, orderOfFiniteElement, deltaT,
							t_end, RUNGE_KUTTA_5STAGE, NONE, is_smooth, semi_lagrangian, false,
							false);

			// Analysis
			// Velocity error (compare Paper by Min and Lee)
			std::stringstream stream1;
			stream1 << "|f-f_ref|_sup; N=" << N << "; p="
					<< orderOfFiniteElement;
			result.quantity.push_back(stream1.str());
			double expected = std::log10(
					400 * std::pow(0.3 * deltaX, orderOfFiniteElement + 1));
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

} /* ConvergenceTestSemiLagrangianAdvection */


TestResult ConvergenceTestSemiLagrangianAdvectionNonsmooth() {

	TestResult result;
	result.id = 13;
	result.name =
			"Convergence Test: SEDG linear advection (non-smooth problem)";
	result.details =
			"This test runs the advection solver for a non-smooth periodic profile, which is "
					"only twice continuously differentiable."
					"The theoretical convergence rate is O(h^min( p+1 , 2 ) * p^(-1.5)) , which is tested.";
	result.time = clock();
	// smooth: 100*(0.3*dx)**(p+1)
	// nonsmooth: (0.5*p)**(-1.5)*4*dx**2

	bool is_smooth = false;
	bool semi_lagrangian = true;
	for (size_t N = 2; N <= 3; N++) {
		for (size_t orderOfFiniteElement = 4; orderOfFiniteElement <= 8;
				orderOfFiniteElement += 4) {

			double deltaX = 1. / (pow(2, N));
			double deltaT = 0.4 * pow(0.5, N)
					/ ((orderOfFiniteElement + 1) * (orderOfFiniteElement + 1));
			double t_end = 0.1;
			if (t_end / deltaT <= 5) {
				continue;
			}
			AdvectionBenchmark::AdvectionResult advectionResult =
					AdvectionBenchmark::oneTest(N, orderOfFiniteElement, deltaT,
							t_end, RUNGE_KUTTA_5STAGE, NONE, is_smooth, semi_lagrangian,
							false);

			// Analysis
			std::stringstream stream1;
			stream1 << "log10(|f-f_ref|_sup); N=" << N << "; p="
					<< orderOfFiniteElement;
			result.quantity.push_back(stream1.str());
			double expected = std::log10(
					0.5 * std::pow(0.5 * orderOfFiniteElement, -1.5) * deltaX
							* deltaX);
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

}
/* ConvergenceSemiLagrangianAdvectionNonsmooth */



} /* namespace IntegrationTests */

