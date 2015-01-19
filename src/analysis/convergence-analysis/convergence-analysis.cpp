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

#ifdef WITH_TRILINOS
	int a = 0;
	char ** b;
	static	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(a,b);
#endif

	/////////////////////////////////////
	// parse command line parameters ////
	/////////////////////////////////////
	bool WALL = false;
	bool DIFFUSIVE = false;

	// default values
	size_t P_MIN = 2;
	size_t P_MAX = 12;
	size_t N_MIN = 2;
	size_t N_MAX = 8;
	double MAX_TIME = 20000; //

	try {
		assert(argc >= 3);
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
		cout
				<< "Optional arguments are (here with default values) -p_min=2 , -p_max=12, -N_min=2, -N_max=8, -max_time=3600"
				<< endl;
		cout
				<< "They define the max and min order of finite element, max and min refinement level and max real time per simulation (in sec)."
				<< endl;
		cout << e.what() << endl;
		return -1;
	}

	// Optional arguments
	size_t max_args = argc;
	for (size_t i = 3; i < max_args; i++) {
		try {
			char delim = '=';
			std::stringstream ss(argv[i]);
			std::string item;
			std::vector<std::string> elems;
			while (std::getline(ss, item, delim)) {
				elems.push_back(item);
			}
			if ("-p_min" == elems[0]) {
				P_MIN = atoi(elems[1].c_str());
			} else if ("-p_max" == elems[0]) {
				P_MAX = atoi(elems[1].c_str());
			} else if ("-N_min" == elems[0]) {
				N_MIN = atoi(elems[1].c_str());
			} else if ("-N_max" == elems[0]) {
				N_MAX = atoi(elems[1].c_str());
			} else if ("-max_time" == elems[0]) {
				MAX_TIME = atof(elems[1].c_str());
			} else {
				cout << "---------------------------" << endl;
				cout << "No option " << elems[0] << endl;
				cout << "---------------------------" << endl;
				throw ("Not all cmd line arguments valid.");
			}
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
			cout
					<< "Optional arguments are (here with default values) -p_min=2 , -p_max=12, -N_min=2, -N_max=8, -max_time=3600"
					<< endl;
			cout
					<< "They define the max and min order of finite element, max and min refinement level and max real time per simulation (in sec)."
					<< endl;
			cout << e.what() << endl;
			return -1;
		}
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
		tmax = L/U*40/(sqrt(3)*20); //40.0 -> dimensionless;

	} else {
		//////////////////////////
		// Taylor-Green vortex ///
		//////////////////////////
		Re =2 * 4 * atan(1); // 2pi
		L = 2 * 4 * atan(1); // 2pi
		U = 1.0; // must not be changed here!
		viscosity = U * L / Re;
		tmax = 0.01 / viscosity; // scale with viscosity to get same behavior for every Re

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
	filename << getenv("NATRIUM_HOME") << "/convergence-analysis/convergence-" << bench_str << "-"
			<< scal_str << "-" << 1900 + ltm->tm_year << "-" << ltm->tm_mon +1
			<< "-" << ltm->tm_mday << "_" << ltm->tm_hour << "-" << ltm->tm_min
			<< ".txt";
	std::ofstream orderFile(filename.str().c_str());
	orderFile
			<< "# refinement   p      dx    #dofs    dt   #steps   tmax    scaling    Ma    tau    max |u_analytic|  ||u_analytic||_2   max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2       init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

	////////////////////////////
	// convergence analysis ////
	////////////////////////////
	cout << "Starting benchmarking for p in [" << P_MIN << ", " << P_MAX << "]; N in [" << N_MIN << ", " << N_MAX << "]; max_time = " << MAX_TIME << "..." << endl;
	for (size_t refinementLevel = N_MIN; refinementLevel <= N_MAX;
			refinementLevel++) {

		for (size_t orderOfFiniteElement = P_MIN; orderOfFiniteElement <= P_MAX;
				orderOfFiniteElement += 1) {
			cout << "N = " << refinementLevel << "; p = "
					<< orderOfFiniteElement << " ... " << endl;

			// make benchmark problem
			if (WALL) {
				shared_ptr<CouetteFlow2D> couette2D =
						make_shared<CouetteFlow2D>(viscosity, U,
								refinementLevel, L, 0.0);
				benchmark = couette2D;
				dx = CFDSolverUtilities::getMinimumVertexDistance<2>(
						*couette2D->getTriangulation());
			} else {
				shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
						TaylorGreenVortex2D>(viscosity, refinementLevel);
				benchmark = tgVortex;
				dx = CFDSolverUtilities::getMinimumVertexDistance<2>(
						*tgVortex->getTriangulation());
			}

			// determine time step size
			if (DIFFUSIVE) {
				double tmp_Ma = Ma * sqrt(dx) / L / orderOfFiniteElement;
				scaling = sqrt(3) * 1 / tmp_Ma;
			}
			double CFL = 0.4;
			dt = CFDSolverUtilities::calculateTimestep<2>(
					*benchmark->getTriangulation(), orderOfFiniteElement,
					D2Q9IncompressibleModel(scaling), CFL);

			// avoid too expensive runs
			// individual jobs should take < 1h
			// estimated runtime /sec : 1e-5 (L / dx * tmax / dt)**(1.5)
			if (1e-5 * pow(L / dx * tmax / dt, 1.5) > MAX_TIME)
				continue;

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
			if ( int(tmax / dt) == 0) {
				continue;
			}
			configuration->setNumberOfTimeSteps(tmax / dt);
			configuration->setInitializationScheme(ITERATIVE);
			configuration->setIterativeInitializationNumberOfIterations(1000);
			configuration->setIterativeInitializationResidual(1e-15);

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
						<< " " << dx << " " << solver.getNumberOfDoFs() << " "
						<< dt << " " << solver.getIteration() << " "
						<< solver.getTime() << " " << scaling << " "
						<< sqrt(3) * U / scaling << " " << solver.getTau()
						<< " " << solver.getErrorStats()->getMaxUAnalytic()
						<< " " << solver.getErrorStats()->getL2UAnalytic()
						<< " " << solver.getErrorStats()->getMaxVelocityError()
						<< " " << solver.getErrorStats()->getMaxDensityError()
						<< " " << solver.getErrorStats()->getL2VelocityError()
						<< " " << solver.getErrorStats()->getL2DensityError()
						<< " " << time1 << " " << time2 << " "
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
