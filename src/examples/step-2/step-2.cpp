/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "CouetteFlow2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-2..." << endl;

	// set Reynolds and Mach number
	const double Re = 10;
	const double Ma = 0.1;

	// set spatial discretization
	size_t refinementLevel = 4;
	size_t orderOfFiniteElement = 2;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1/sqrt(3)*Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)

	// in order to start from a continuous solution, do not start at t=0
	const double startTime = 0.0;

	// set small time step size
	const double timeStepSize = 0.0001;

	cout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;
	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-2";
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setNumberOfTimeSteps(5./timeStepSize);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(timeStepSize);
	configuration->setCommandLineVerbosity(7);
	//configuration->setDistributionInitType(Iterative);

	shared_ptr<CouetteFlow2D> couetteFlow = make_shared<CouetteFlow2D>(
			viscosity, U, refinementLevel, 1.0, startTime);
	shared_ptr<Benchmark<2> > couetteProblem = couetteFlow;
	BenchmarkCFDSolver<2> solver(configuration, couetteProblem);

	solver.run();

	cout << "NATriuM step-2 terminated." << endl;

	return 0;
}
