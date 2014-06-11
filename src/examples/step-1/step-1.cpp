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

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-1 ..." << endl;

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

	const double refinementLevel = 5;
		cout << "refinement Level = " << refinementLevel << endl;
//		for (size_t orderOfFiniteElement = 2; orderOfFiniteElement < 7;
//				orderOfFiniteElement++) {
//			cout << "FE order = " << orderOfFiniteElement << endl;
		// the scaling has to be orders of magnitude greater than the boundary velocity

		// calculate distance between quadrature nodes
		double dx = 2 * 3.1415926
				/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));
		// chose dt so that courant (advection) = 1 for the diagonal directions
		double dt = dx / (scaling * sqrt(2));

		cout << "dt = " << dt << " ...";

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/step-1";
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
		configuration->setCommandLineVerbosity(BASIC);
		configuration->setTimeStepSize(dt);
		if (dt > 0.1) {
			cout << "Timestep too big." << endl;
		}
//configuration->setNumberOfTimeSteps(1.0 / dt);
		configuration->setNumberOfTimeSteps(1);

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
		} catch (std::exception& e) {
			cout << " Error" << endl;
		}

	cout << "step-1 terminated." << endl;

	return 0;
}
