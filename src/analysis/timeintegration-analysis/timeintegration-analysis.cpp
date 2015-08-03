/**
 * @file convergence-analysis-basic.cpp
 * @short The convergence of the NATriuM solver is analyzed by application to the Taylor-Green vortex in 2D (only periodic walls).
 * This script uses a linear scaling (= constant Mach number). Thus, there is a general compressibily error, which destroys the
 * convergence for the finer meshes. To analyze the results, move the table_order.txt and table_results.txt files to NATriuM/src/analysis/convergence_analysis_basic/
 * and execute the gnuplot scripts.
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */
/**
 * @file convergence-analysis-p.cpp
 * @short The convergence of the NATriuM solver is analyzed by application to the Taylor-Green vortex in 2D (only periodic walls).
 * This script uses a linear scaling (= constant Mach number). Thus, there is a general compressibily error, which destroys the
 * convergence for the finer meshes. To analyze the results, move the table_order.txt and table_results.txt
 * files to NATriuM/src/analysis/convergence_analysis_basic/
 * and execute the gnuplot scripts.
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <string>

#include "deal.II/numerics/data_out.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/CouetteFlow2D.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// if this define statement is enabled: only the initialization time is regarded
//#define MEASURE_ONLY_INIT_TIME
// #define ONLY_PERIODIC

// Main function
int main() {

	cout
			<< "Starting NATriuM time integrator analysis with moving wall boundaries..."
			<< endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

#ifdef ONLY_PERIODIC
	// specify Reynolds number
	const double Re = 8 * atan(1); // = 2 pi
	// specify Mach number
	const double Ma = 0.1;
#else
	const double Re = 2000;
	const double Ma = 0.05;
#endif

	// Problem description
	// -------------------
#ifdef ONLY_PERIODIC
	// length of quadratic domain
	const double L =  8 * atan(1); // = 2 pi
	const double U = 1;
	const double tmax = 1;
#else
	const double L = 1;
	// velocity of top plate
	const double U = 1 / sqrt(3) * Ma;//5.773502691896258e-02/3.1415926;//0.02;
	const double tmax = 10;
#endif
	// scaling of particle velocities
	double scaling = sqrt(3) * U / Ma;
	// Viscosity
	const double viscosity = U * L / Re;
	// starting time
	//const double t0 = 30.0;
	const double t0 = 1.0; 	// analytic solution won't converge for t0 = 0.0 and adaptive timesteps

	size_t refinementLevel = 3;
	size_t orderOfFiniteElement = 5;

#ifndef ONLY_PERIODIC
	// make problem object
	shared_ptr<CouetteFlow2D> couette2D = make_shared<CouetteFlow2D>(viscosity,
			U, refinementLevel, L, t0, false);
	shared_ptr<Benchmark<2> > benchmark = couette2D;
#else
	// make problem object
	shared_ptr<TaylorGreenVortex2D> tgv2D = make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel, 1. / sqrt(3.) / Ma);
	shared_ptr<Benchmark<2> > benchmark = tgv2D;
#endif

	// prepare time table file
	// the output is written to the standard output directory (e.g. NATriuM/results or similar)

	for (int integrator = 1; integrator < 10; integrator++) {
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		std::string timeintegrator = "test";

		switch (integrator) {
		case 1:
			configuration->setTimeIntegrator(RUNGE_KUTTA_5STAGE);
			timeintegrator = "RUNGE_KUTTA_5STAGE";
			break;

		case 2:
			configuration->setTimeIntegrator(THETA_METHOD);
			timeintegrator = "THETA_METHOD";
			break;

		case 3:
			configuration->setTimeIntegrator(EXPONENTIAL);
			timeintegrator = "EXPONENTIAL";
			break;

		case 4:
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(FORWARD_EULER);
			timeintegrator = "FORWARD_EULER";
			break;

		case 5:
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(RK_THIRD_ORDER);
			timeintegrator = "RK_THIRD_ORDER";
			break;

		case 6: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(RK_CLASSIC_FOURTH_ORDER);
			timeintegrator = "RK_CLASSIC_FOURTH_ORDER";
			break;
		}
		case 7: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(BACKWARD_EULER);
			timeintegrator = "BACKWARD_EULER";
			break;
		}
		case 8: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(IMPLICIT_MIDPOINT);
			timeintegrator = "IMPLICIT_MIDPOINT";
			break;
		}
		case 9: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(CRANK_NICOLSON);
			timeintegrator = "CRANK_NICOLSON";
			break;
		}
		case 10: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(SDIRK_TWO_STAGES);
			timeintegrator = "SDIRK_TWO_STAGES";
			break;
		}
		case 11: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(HEUN_EULER);
			timeintegrator = "HEUN_EULER";
			break;
		}
		case 12: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(BOGACKI_SHAMPINE);
			timeintegrator = "BOGACKI_SHAMPINE";
			break;
		}
		case 13: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(DOPRI);
			timeintegrator = "DOPRI";
			break;
		}
		case 14: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(FEHLBERG);
			timeintegrator = "FEHLBERG";
			break;
		}
		case 15: {
			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(CASH_KARP);
			timeintegrator = "CASH_KARP";
			break;
		}
		}

		std::stringstream filename;
		filename << getenv("NATRIUM_HOME") << "/timeintegration-analysis/"
				<< timeintegrator.c_str() << "_table.txt";
		std::ofstream timeFile(filename.str().c_str());
		timeFile
				<< "# order of FE   dt        init time (sec)             loop time (sec)         time for one iteration (sec)"
				<< endl;

		// prepare error table file
		std::stringstream filename2;
		filename2 << getenv("NATRIUM_HOME") << "/timeintegration-analysis/"
				<< timeintegrator.c_str() << "_error.txt";
		std::ofstream orderFile(filename2.str().c_str());
		orderFile << "# visc = " << viscosity << "; Ma = " << Ma << endl;
		orderFile
				<< "#  dt  i      CFL         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   "
						"||error_rho||_2	runtime" << endl;

		cout << "Time integrator: " << timeintegrator.c_str() << endl;

		for (double CFL = 0.2; CFL <= 102.4; CFL *= 2.) {

			double dt = CFDSolverUtilities::calculateTimestep<2>(
					*benchmark->getMesh(), orderOfFiniteElement,
					D2Q9(scaling), CFL);

			cout << "CFL = " << CFL << endl;
			if (tmax/ dt < 2) {
				cout << "time step too big." << endl;
				continue;
			}

			// time measurement variables
			double time1, time2, timestart;

			// setup configuration
			std::stringstream dirName;
			dirName << getenv("NATRIUM_HOME") << "/timeintegration-analysis/"
					<< timeintegrator << "_" << CFL << "_"
					<< dt;

			//configuration->setSwitchOutputOff(true);
			configuration->setOutputDirectory(dirName.str());
			configuration->setRestartAtLastCheckpoint(false);
			configuration->setUserInteraction(false);
			configuration->setOutputTableInterval(10000000);

			//configuration->setOutputCheckpointInterval(1000);
			configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
			configuration->setStencilScaling(scaling);
			configuration->setCommandLineVerbosity(BASIC);
			configuration->setTimeStepSize(dt);
			configuration->setSimulationEndTime(tmax);

#ifdef MEASURE_ONLY_INIT_TIME
			configuration->setNumberOfTimeSteps(1);
#endif

			// solver
			timestart = clock();
			BenchmarkCFDSolver<2> solver(configuration, benchmark);
			time1 = clock() - timestart;

			try {
				solver.run();
				time2 = clock() - time1 - timestart;
				time1 /= CLOCKS_PER_SEC;
				time2 /= CLOCKS_PER_SEC;
				solver.getErrorStats()->update();

				/*if (solver.getErrorStats()->getMaxVelocityError()
						> solver.getErrorStats()->getMaxUAnalytic())
					break;*/

				cout << " OK ... Init: " << time1 << " sec; Run: " << time2
						<< " sec." << " Max Error U_analytic: "
						<< solver.getErrorStats()->getMaxVelocityError()
						<< endl;
				// put out runtime
				timeFile << orderOfFiniteElement << "         " << dt
						<< "      " << time1 << "     " << time2 << "        "
						<< time2 / configuration->getNumberOfTimeSteps()
						<< endl;
				// put out final errors
				solver.getErrorStats()->update();
				orderFile << dt << " " << solver.getIteration() << " "
						<< CFL << " "
						<< solver.getErrorStats()->getMaxUAnalytic() << " "
						<< solver.getErrorStats()->getMaxVelocityError() << " "
						<< solver.getErrorStats()->getMaxDensityError() << " "
						<< solver.getErrorStats()->getL2VelocityError() << " "
						<< solver.getErrorStats()->getL2DensityError() << " "
						<< time2 << endl;
			} catch (std::exception& e) {
				cout << " Error: " << e.what() << endl;
			}

		} /* for time step*/
	} /* for integrator */
	cout << "Time integration analysis terminated." << endl;

	return 0;
}
