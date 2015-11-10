/**
 * @file step-1.cpp
 * @short Taylor-Green vortex in 2D (only periodic walls)
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	MPIGuard::getInstance();

	pout << "Starting NATriuM step-1 ..." << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	const double L = 8 * atan(1);
	const double U = 1;
	const double Ma = 0.2;
	const double refinementLevel = 1;
	double scaling = sqrt(3) * 1 / Ma;
	const bool iterativeInitialization = false;

	for (double Re = 0.1*8*atan(1); Re <= 10.0; Re*=10) {

		// Re = viscosity/(2*pi)
		const double viscosity = L * U / Re;
		// zunaechst: fixed order of FE
		for (double orderOfFiniteElement = 2; orderOfFiniteElement <= 10;
				orderOfFiniteElement += 2) {

			// make problem object
			shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
					TaylorGreenVortex2D>(viscosity, refinementLevel);
			shared_ptr<Benchmark<2> > taylorGreen = tgVortex;

			double dx = CFDSolverUtilities::getMinimumVertexDistance<2>(
									*tgVortex->getMesh());
			double tmp_Ma = Ma * sqrt(dx) / L / orderOfFiniteElement;
			scaling = sqrt(3) * 1 / tmp_Ma;

			double dt = CFDSolverUtilities::calculateTimestep<2>(
					*tgVortex->getMesh(), orderOfFiniteElement,
					D2Q9IncompressibleModel(scaling), 0.4);
			pout << "p = " << orderOfFiniteElement << "; dt = " << dt << " ..."
					<< endl;

			// time measurement variables
			double time1, time2, timestart;

			// setup configuration
			std::stringstream dirName;
			dirName << getenv("NATRIUM_HOME") << "/step-1-variousRe/";
			if (iterativeInitialization) {
				dirName << "iter_init-";
			} else {
				dirName << "eq_init-";
			}
			dirName << "Re" << Re << "-p" << orderOfFiniteElement;
			shared_ptr<SolverConfiguration> configuration = make_shared<
					SolverConfiguration>();
			//configuration->setSwitchOutputOff(true);
			configuration->setOutputDirectory(dirName.str());
			configuration->setRestartAtLastCheckpoint(false);
			configuration->setUserInteraction(false);
			configuration->setOutputTableInterval(1);
			configuration->setOutputCheckpointInterval(1e9);
			configuration->setOutputSolutionInterval(1e9);
			configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
			configuration->setStencilScaling(scaling);
			configuration->setCommandLineVerbosity(WARNING);
			configuration->setTimeStepSize(dt);

			if (iterativeInitialization) {
				configuration->setInitializationScheme(ITERATIVE);
				configuration->setIterativeInitializationNumberOfIterations(
						1000);
				configuration->setIterativeInitializationResidual(1e-15);
			}

			if (dt > 0.1) {
				pout << "Timestep too big." << endl;
			}

			configuration->setNumberOfTimeSteps(2.0 / dt);

			timestart = clock();
			BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
			time1 = clock() - timestart;

			try {
				solver.run();
				time2 = clock() - time1 - timestart;
				time1 /= CLOCKS_PER_SEC;
				time2 /= CLOCKS_PER_SEC;
				pout << " OK ... Init: " << time1 << " sec; Run: " << time2
						<< " sec." << endl;
			} catch (std::exception& e) {
				pout << " Error" << endl;
			}

			LOG(BASIC) << "NATriuM run for Re " << Re << " complete." << endl;
		}
	}
	pout << "step-1-variousRe terminated." << endl;

	return 0;
}
