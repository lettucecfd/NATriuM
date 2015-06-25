/**
 * @file convergence-analysis-basic.cpp
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

#include <boost/filesystem.hpp>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/Logging.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// if this define statement is enabled: only the initialization time is regarded
// #define MEASURE_ONLY_INIT_TIME

// Main function
int main(int argc, char* argv[]) {

	cout << "Starting NATriuM convergence analysis (Ma-dependency)..." << endl;

	// PARSE INPUT
	cout
			<< "To use incompressible scheme, pass 1st cmd line parameter 1 (standard), 2 (trafo) or 3 (incompressible)."
			<< endl;
	cout
			<< "To specify initialization, pass 2nd cmd line parameter 1 (rho=1), 2 (iterative) or 3 (analytic)."
			<< endl;
	cout
			<< "To specify a CFL number other than 0.4, pass CFL number as 3rd cmd line parameter."
			<< endl;
	int BGK_SCHEME = 1;
	int INIT_SCHEME = 1;
	double CFL = 0.4;
	if (argc > 0) {
		if (std::atoi(argv[1]) == 1) {
			cout << "Standard BGK" << endl;
			BGK_SCHEME = 1;
		} else if (std::atoi(argv[1]) == 2) {
			cout << "Transformed BGK" << endl;
			BGK_SCHEME = 2;
		} else if (std::atoi(argv[1]) == 3) {
			cout << "Incompressible scheme by He and Luo" << endl;
			BGK_SCHEME = 3;
		} else {
			cout << "Did not understand collision scheme." << endl;
		}
	}
	if (argc > 1) {
		if (std::atoi(argv[2]) == 1) {
			cout << "Init with rho = 1" << endl;
			INIT_SCHEME = 1;
		} else if (std::atoi(argv[2]) == 2) {
			cout << "Init with iterative scheme" << endl;
			INIT_SCHEME = 2;
		} else if (std::atoi(argv[2]) == 3) {
			cout << "Init with analytic pressure" << endl;
			INIT_SCHEME = 3;
		} else {
			cout << "Did not understand init scheme." << endl;
		}
	}
	if (argc > 2) {
		CFL = std::atof(argv[3]);
	}

/////////////////////////////////////////////////
// set parameters, set up configuration object
//////////////////////////////////////////////////

// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/convergence-analysis-Ma/";
	if (1 == BGK_SCHEME) {
		dirName << "std";
	} else if (2 == BGK_SCHEME) {
		dirName << "trans";
	} else if (3 == BGK_SCHEME) {
		dirName << "inc";
	}
	if (1 == INIT_SCHEME) {
		dirName << "_rho1/";
	} else if (2 == INIT_SCHEME) {
		dirName << "_iter/";
	} else if (3 == INIT_SCHEME) {
		dirName << "_ana/";
	}
	boost::filesystem::create_directory(dirName.str().c_str());

	const double viscosity = 1.;

// zunaechst: fixed order of FE
	const size_t refinementLevel = 4;

// prepare time table file
// the output is written to the standard output directory (e.g. NATriuM/results or similar)
	std::stringstream filename;
	filename << dirName.str() << "table_runtime.txt";
	std::ofstream timeFile(filename.str().c_str());
	timeFile
			<< "#refinement Level    FE order     dt        init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

// prepare error table file
	std::stringstream filename2;
	filename2 << dirName.str() << "table_order.txt";
	std::ofstream orderFile(filename2.str().c_str());
//orderFile << "# visc = " << viscosity << "; Ma = " << Ma << endl;
	orderFile
			<< "#  refinementlevel  FE order    Ma    dt    tau    i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
			<< endl;

	for (size_t orderOfFiniteElement = 1; orderOfFiniteElement < 5;
			orderOfFiniteElement++) {
		cout << "refinement Level = " << refinementLevel << endl;

		for (double Ma = 0.3; Ma > 5e-4; Ma /= 2) {
			cout << "Ma = " << Ma << endl;

			double scaling = sqrt(3) * 1 / Ma;

			shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
					TaylorGreenVortex2D>(viscosity, refinementLevel,
					scaling / sqrt(3), (INIT_SCHEME == 3));
			shared_ptr<Benchmark<2> > taylorGreen = tgVortex;

			double dt = CFDSolverUtilities::calculateTimestep<2>(
					*tgVortex->getTriangulation(), orderOfFiniteElement,
					D2Q9(scaling), CFL);

			cout << "dt = " << dt << " ...";

			// time measurement variables
			double time1, time2, timestart;

			// setup configuration
			std::stringstream dirName2;
			dirName2 << dirName.str() << Ma << "-ref" << refinementLevel << "-p"
					<< orderOfFiniteElement;
			shared_ptr<SolverConfiguration> configuration = make_shared<
					SolverConfiguration>();
			//configuration->setSwitchOutputOff(true);
			configuration->setOutputDirectory(dirName2.str());
			configuration->setRestartAtLastCheckpoint(false);
			configuration->setUserInteraction(false);
			configuration->setOutputTableInterval(100);
			configuration->setOutputSolutionInterval(100000);
			configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
			configuration->setStencilScaling(scaling);
			configuration->setCommandLineVerbosity(0);
			if (2 == INIT_SCHEME) {
				configuration->setInitializationScheme(ITERATIVE);
				configuration->setIterativeInitializationNumberOfIterations(
						1e5);
				configuration->setIterativeInitializationResidual(1e-10);

			}
			if (1 == BGK_SCHEME) {
				configuration->setCollisionScheme(BGK_STANDARD);
			}
			if (2 == BGK_SCHEME) {
				configuration->setCollisionScheme(BGK_STANDARD_TRANSFORMED);
			} else if (3 == BGK_SCHEME) {
				cout << "Not yet implemented" << endl;
				continue;
			}
			configuration->setTimeStepSize(dt);
			if (dt > 0.1) {
				cout << "Timestep too big." << endl;
				continue;

			}
			configuration->setNumberOfTimeSteps(1.0 / dt);

#ifdef MEASURE_ONLY_INIT_TIME
			configuration->setNumberOfTimeSteps(1);
#endif

			// make problem and solver objects; measure time
			timestart = clock();
			BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
			time1 = clock() - timestart;

			try {
				solver.run();
				time2 = clock() - time1 - timestart;
				time1 /= CLOCKS_PER_SEC;
				time2 /= CLOCKS_PER_SEC;
				cout << " OK ... Init: " << time1 << " sec; Run: " << time2
						<< " sec." << endl;
				// put out runtime
				timeFile << refinementLevel << "     " << orderOfFiniteElement
						<< "         " << dt << "      " << time1 << "     "
						<< time2 << "        "
						<< time2 / configuration->getNumberOfTimeSteps()
						<< endl;
				// put out final errors
				solver.getErrorStats()->update();
				orderFile << refinementLevel << "     " << orderOfFiniteElement
						<< " " << Ma << "  " << dt << " " << solver.getTau()
						<< " " << solver.getIteration() << " "
						<< solver.getTime() << " "
						<< solver.getErrorStats()->getMaxUAnalytic() << " "
						<< solver.getErrorStats()->getMaxVelocityError() << " "
						<< solver.getErrorStats()->getMaxDensityError() << " "
						<< solver.getErrorStats()->getL2VelocityError() << " "
						<< solver.getErrorStats()->getL2DensityError() << endl;
			} catch (std::exception& e) {
				cout << " Error" << endl;
			}
		} /* for Ma */

		orderFile << endl;
		timeFile << endl;

	} /* for refinement level */

	cout << "Convergence analysis (Ma) terminated." << endl;

	return 0;
}
