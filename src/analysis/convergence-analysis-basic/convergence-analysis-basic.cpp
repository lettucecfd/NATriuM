/**
 * @file convergence-analysis-basic.cpp
 * @short Taylor-Green vortex in 2D (only periodic walls)
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "../../examples/step-1/TaylorGreenVortex2D.h"

using namespace natrium;

// if this define statement is enabled: only the initialization time is regarded
//#define MEASURE_ONLY_INIT_TIME


// Main function
int main() {

	cout << "Starting NATriuM convergence analysis (linear scaling)..." << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	const double viscosity = 1;
	// C-E-approach: constant stencil scaling
	// specify Mach number
	const double Ma = 0.05;
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = 2;

	// chose scaling so that the right Ma-number is achieved
	double scaling = sqrt(3) * 1 / Ma;

	// prepare table file
	std::stringstream filename;
	filename << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-basic/runtime.txt";
	std::ofstream timeFile(filename.str().c_str());
	timeFile
			<< "#refinement Level     dt        init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

	for (size_t refinementLevel = 2; refinementLevel < 9; refinementLevel++) {
		cout << "refinement Level = " << refinementLevel << endl;

		double dx = 2 * 3.1415926
				/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));
		// chose dt so that courant (advection) = 1 for the diagonal directions
		double dt = dx / (scaling * sqrt(2));

		cout << "dt = " << dt << " ...";

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << "../results/convergence-analysis-basic/"
				<< orderOfFiniteElement << "_" << refinementLevel;
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(10);
		configuration->setOutputCheckpointInterval(1000);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(0);
		configuration->setTimeStepSize(dt);
		if (dt > 0.1) {
			cout << "Timestep too big." << endl;
			continue;

		}
		configuration->setNumberOfTimeSteps(1.0 / dt);

#ifdef MEASURE_ONLY_INIT_TIME
		configuration->setNumberOfTimeSteps(1);
#endif

		// make problem and solver objects
		shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel);
		shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
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
			timeFile << refinementLevel << "         " << dt << "      "
					<< time1 << "     " << time2 << "        " << time2/configuration->getNumberOfTimeSteps() << endl;
		} catch (std::exception& e) {
			cout << " Error" << endl;
		}
	}
	cout << "Convergence analysis (basic) terminated." << endl;

	return 0;
}
