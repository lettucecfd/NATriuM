/**
 * @file performance-analysis.cpp
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

#include "natrium/stencils/D2Q9.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	MPIGuard::getInstance();

	pout << "Starting NATriuM performance analysis..." << endl;

	/////////////////////////////////////
	// parse command line parameters ////
	/////////////////////////////////////
	bool WALL = false;
	bool DIFFUSIVE = false;

	// default values
	size_t P_MIN = 2;
	size_t P_MAX = 10;
	size_t N_MIN = 2;
	size_t N_MAX = 7;

	// Optional arguments
	size_t max_args = argc;
	for (size_t i = 1; i < max_args; i++) {
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
			} else {
				pout << "---------------------------" << endl;
				pout << "No option " << elems[0] << endl;
				pout << "---------------------------" << endl;
				throw("Not all cmd line arguments valid.");
			}
		} catch (std::exception& e) {
			pout
					<< "Optional arguments are (here with default values) -p_min=2 , -p_max=12, -N_min=2, -N_max=8"
					<< endl;
			pout
					<< "They define the max and min order of finite element + max and min refinement level."
					<< endl;
			pout << e.what() << endl;
			return -1;
		}
	}

	///////////////////////////////
	// set flow characteristics ///
	///////////////////////////////

	// Benchmark problem
	boost::shared_ptr<ProblemDescription<2> > benchmark;
	// dt
	double dt;
	// number of time steps
	size_t number_of_steps = 20;

	//////////////////////////
	// Taylor-Green vortex ///
	//////////////////////////
	double Re = 2 * 4 * atan(1); // 2pi
	double L = 2 * 4 * atan(1); // 2pi
	double U = 1.0; // must not be changed here!
	double viscosity = U * L / Re;
	double Ma = 0.05; // constant Ma number
	double scaling = sqrt(3) * 1 / Ma;

	///////////////////////////
	// prepare table file  ////
	///////////////////////////
	time_t now = time(0);
	tm *ltm = localtime(&now);
	std::stringstream filename_ntrm;
	filename_ntrm << getenv("NATRIUM_HOME")
			<< "/performance-analysis/performance-NTrM-" << getenv("HOSTNAME")
			<< "-" << 1900 + ltm->tm_year << "-" << ltm->tm_mon + 1 << "-"
			<< ltm->tm_mday << "_" << ltm->tm_hour << "-" << ltm->tm_min
			<< ".txt";
	std::ofstream outfile_ntrm(filename_ntrm.str().c_str());
	outfile_ntrm
			<< "# refinement   p    #dofs   #steps    init time (sec)   loop time (sec)   time for one iteration (sec)"
			<< endl;

	////////////////////////////////////
	// performance analysis natrium ////
	////////////////////////////////////
	pout << "Starting Benchmarking for p in [" << P_MIN << ", " << P_MAX
			<< "]; N in [" << N_MIN << ", " << N_MAX << "]" << endl;
	for (size_t orderOfFiniteElement = P_MIN; orderOfFiniteElement <= P_MAX;
			orderOfFiniteElement += 1) {

		for (size_t refinementLevel = N_MIN; refinementLevel <= N_MAX;
				refinementLevel++) {
			pout << "N = " << refinementLevel << "; p = "
					<< orderOfFiniteElement << " ... " << endl;
			// make benchmark problem

			boost::shared_ptr<TaylorGreenVortex2D> tgVortex = boost::make_shared<
					TaylorGreenVortex2D>(viscosity, refinementLevel);
			benchmark = tgVortex;

			double CFL = 0.4;

			/////////////////////////////
			// run benchmark problem ////
			/////////////////////////////

			// time measurement variables
			double time1, time2, timestart;

			// setup configuration
			std::stringstream dirName;
			boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
					SolverConfiguration>();
			configuration->setSwitchOutputOff(true);
			configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
			configuration->setStencilScaling(scaling);
			configuration->setCommandLineVerbosity(WARNING);
			configuration->setCFL(CFL);
			configuration->setNumberOfTimeSteps(20);
			timestart = clock();
			CFDSolver<2> solver(configuration, benchmark);
			time1 = clock() - timestart;

			size_t n_steps = 20;
			try {
				for (size_t i = 0; i < n_steps; i++) {
					solver.stream();
					solver.collide();
				}
				time2 = clock() - time1 - timestart;
				time1 /= CLOCKS_PER_SEC;
				time2 /= CLOCKS_PER_SEC;
				pout << " OK ... Init: " << time1 << " sec; Run: " << time2
						<< " sec." << endl;

				// put out errors and times
				outfile_ntrm << refinementLevel << " " << orderOfFiniteElement
						<< " " << solver.getNumberOfDoFs() << " "
						<< n_steps << " " << " " << time1 << " "
						<< time2 << " "
						<< time2 / configuration->getNumberOfTimeSteps()
						<< endl;
			} catch (std::exception& e) {
				pout << " Error: " << e.what() << endl;
			}

		} /* for order FE */
		outfile_ntrm << endl;
	} /* for refinement level */

	pout << "Performance analysis terminated." << endl;
	return 0;
}
