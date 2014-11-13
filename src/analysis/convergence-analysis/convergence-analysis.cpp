/**
 * @file convergence-analysis.cpp
 * @short The convergence of the NATriuM solver is analyzed by application to the Taylor-Green vortex in 2D (only periodic walls).
 * This script uses a linear scaling (= constant Mach number). Thus, there is a general compressibily error, which destroys the
 * convergence for the finer meshes. To analyze the results, move the table_order.txt and table_results.txt files to NATriuM/src/analysis/convergence_analysis_basic/
 * and execute the gnuplot scripts.
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <exception>
#include <ctime>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"
#include "utilities/CFDSolverUtilities.h"

#include "../../examples/step-1/TaylorGreenVortex2D.h"
#include "../../examples/step-2/CouetteFlow2D.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	cout << "Starting NATriuM convergence analysis (h + p)..." << endl;

	/////////////////////////////////////
	// parse command line parameters ////
	/////////////////////////////////////
	bool WALL = false;
	bool DIFFUSIVE = false;

	try {
		assert(argc == 3);
		size_t wall = atoi(argv[1]);
		size_t diffusive = atoi(argv[2]);
		assert(wall < 2);
		assert(diffusive < 2);
		WALL = (wall == 1);
		DIFFUSIVE = (diffusive == 1);
	} catch (std::exception& e) {
		cout
				<< "Please call this program with command line arguments( e.g. convergence-analysis-hp 0 0)."
				<< endl;
		cout
				<< "The first number (0 or 1) indicates whether you want to use wall boundaries (i.e. unsteady Couette flow benchmark)."
				<< endl;
		cout
				<< "The second number (0 or 1) indicates whether you want to use the diffusive scaling (i.e. dt ~ dx^2)."
				<< endl;
		cout << e.what() << endl;
		return -1;
	}

	///////////////////////////////
	// set flow characteristics ///
	///////////////////////////////

	// Mach number
	double Ma;
	// Reynolds number
	double Re;
	// characteristic velocity
	double U;
	// characteristic length
	double L;
	// viscosity
	double viscosity;
	// scaling of stencil
	double scaling;
	// Benchmark problem
	shared_ptr<Benchmark<2> > benchmark;
	// dx
	double dx;
	// dt
	double dt;
	// end time
	double tmax;

	if (WALL) {
		//////////////////////////
		// Couette Flow //////////
		//////////////////////////
		Re = 2000;
		U = 1;
		L = 1;
		viscosity = U * L / Re;
		tmax = 40.0;

	} else {
		//////////////////////////
		// Taylor-Green vortex ///
		//////////////////////////
		Re = 2 * 4 * atan(1); // 2pi
		L = 2 * 4 * atan(1); // 2pi
		U = 1.0; // must not be changed here!
		viscosity = U * L / Re;
		tmax = 1.0;

	}

	if (DIFFUSIVE) {
		Ma = 0.2; // Ma number for coarsest discretization
		scaling = sqrt(3) * 1 / Ma;
	} else {
		Ma = 0.05; // constant Ma number
		scaling = sqrt(3) * 1 / Ma;
	}

	///////////////////////////
	// prepare table file  ////
	///////////////////////////

	string bench_str = WALL ? "wall" : "periodic";
	string scal_str = DIFFUSIVE ? "diffusive" : "acoustic";
	time_t now = time(0);
	tm *ltm = localtime(&now);
	std::stringstream filename;
	filename << getenv("NATRIUM_HOME") << "/convergence-" << bench_str << "-"
			<< scal_str << "-" << 1900 + ltm->tm_year << "-" << ltm->tm_mon	<< "-" << ltm->tm_mday << "_" << ltm->tm_hour << "-" << ltm->tm_min << ".txt";
	std::ofstream orderFile(filename.str().c_str());
	orderFile
			<< "# refinement   p      dx    #dofs    dt   #steps   tmax    scaling    Ma    tau    max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2       init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

	////////////////////////////
	// convergence analysis ////
	////////////////////////////
	for (size_t refinementLevel = 2; refinementLevel <= 10; refinementLevel++) {

		for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 14;
				orderOfFiniteElement += 2) {
			cout << "N = " << refinementLevel << "; p = "
					<< orderOfFiniteElement << " ... " << endl;

			// make benchmark problem
			if (WALL) {
				shared_ptr<CouetteFlow2D> couette2D =
						make_shared<CouetteFlow2D>(viscosity, U,
								refinementLevel, L, 0.0);
				benchmark = couette2D;
				dx = CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
						*couette2D->getTriangulation(), orderOfFiniteElement);
			} else {
				shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
						TaylorGreenVortex2D>(viscosity, refinementLevel);
				benchmark = tgVortex;
				dx = CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
						*tgVortex->getTriangulation(), orderOfFiniteElement);
			}

			// determine time step size
			if (DIFFUSIVE) {
				double tmp_Ma = Ma * dx / L;
				scaling = sqrt(3) * 1 / tmp_Ma;
			}
			double CFL = 0.4;
			dt = CFDSolverUtilities::calculateTimestep<2>(
					*benchmark->getTriangulation(), orderOfFiniteElement,
					D2Q9IncompressibleModel(scaling), CFL);

			// avoid too expensive runs
			// individual jobs should take < 1h
			if (L/dx * tmax/dt > 3 * 1e3*1e4) {
				continue;
			}

			/////////////////////////////
			// run benchmark problem ////
			/////////////////////////////

			// time measurement variables
			double time1, time2, timestart;

			// setup configuration
			std::stringstream dirName;
			shared_ptr<SolverConfiguration> configuration = make_shared<
					SolverConfiguration>();
			configuration->setSwitchOutputOff(true);
			configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
			configuration->setStencilScaling(scaling);
			configuration->setCommandLineVerbosity(WARNING);
			configuration->setTimeStepSize(dt);
			configuration->setNumberOfTimeSteps(tmax / dt);

			timestart = clock();
			BenchmarkCFDSolver<2> solver(configuration, benchmark);
			time1 = clock() - timestart;

			try {
				solver.run();
				time2 = clock() - time1 - timestart;
				time1 /= CLOCKS_PER_SEC;
				time2 /= CLOCKS_PER_SEC;
				cout << " OK ... Init: " << time1 << " sec; Run: " << time2
						<< " sec." << endl;

				// put out errors and times
				solver.getErrorStats()->update();
				orderFile << refinementLevel << " " << orderOfFiniteElement
						<< " " << dx << " " << solver.getNumberOfDoFs() << " " << dt << " "
						<< solver.getIteration() << " " << solver.getTime()
						<< " " << scaling << " " << sqrt(3) * U / scaling << " "
						<< solver.getTau() << " "
						<< solver.getErrorStats()->getMaxUAnalytic() << " "
						<< solver.getErrorStats()->getMaxVelocityError() << " "
						<< solver.getErrorStats()->getMaxDensityError() << " "
						<< solver.getErrorStats()->getL2VelocityError() << " "
						<< solver.getErrorStats()->getL2DensityError() << " "
						<< time1 << " " << time2 << " "
						<< time2 / configuration->getNumberOfTimeSteps()
						<< endl;
			} catch (std::exception& e) {
				cout << " Error: " << e.what() << endl;
			}

		} /* for order FE */
		orderFile << endl;
	} /* for refinement level */

	cout << "Convergence analysis terminated." << endl;

	return 0;
}
