/**
 * @file convergence-analysis-kinE.cpp
 * @short Investigate the evolution of the kinetic energy over time
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/collision/BGK.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	MPIGuard::getInstance();

	pout << "Starting Taylor Green vortex to check the kinetic energy decay..." << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	const double viscosity = 1;
	//initial mach number = 0.05
	const double Ma = 0.05;
	// fixed order of FE and refinement level
	const double orderOfFiniteElement = 2;
	const size_t refinementLevel = 3;

	// chose scaling so that Ma is recovered
	double scaling = sqrt(3) * 1 / Ma;

	boost::shared_ptr<TaylorGreenVortex2D> tgVortex = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel);

	// chose dt so that courant (advection) = 1 for the diagonal directions
	double CFL = 0.4;
	double 	dt = CFDSolverUtilities::calculateTimestep<2>(
			*tgVortex->getMesh(), orderOfFiniteElement,
			D2Q9(scaling), CFL);



	// time measurement variables
	//double time1, time2, timestart;

	// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/convergence-analysis-kinE";
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(10);
	configuration->setOutputSolutionInterval(1000);
	configuration->setOutputCheckpointInterval(1000);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(0);
	configuration->setCFL(CFL);
	configuration->setConvergenceThreshold(1e-15);
	if (dt > 0.1) {
		pout << "Timestep too big." << endl;

	}
	configuration->setNumberOfTimeSteps(20.0 / dt);

	// make problem and solver objects
	boost::shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
	//timestart = clock();
	BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
	// get tau to test if the "constant value" is really constant
	const double tau = BGK::calculateRelaxationParameter(viscosity,
			dt, *solver.getStencil());
	pout << "... scaling = " << scaling << " ... tau = " << tau << " ..." << endl;
	//time1 = clock() - timestart;

	try {
		solver.run();
		/*time2 = clock() - time1 - timestart;
		time1 /= CLOCKS_PER_SEC;
		time2 /= CLOCKS_PER_SEC;
		pout << " OK ... Init: " << time1 << " sec; Run: " << time2 << " sec."
				<< endl;*/
	} catch (std::exception& e) {
		pout << " Error" << endl;
	}

pout << "Kinetic energy decay test terminated." << endl;

return 0;
}
