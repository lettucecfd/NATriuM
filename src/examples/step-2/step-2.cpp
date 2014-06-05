/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "CouetteFlow2D.h"

//#define PRINT_SYSTEM_VECTOR

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-2..." << endl;

	// set parameters, set up configuration object
	size_t refinementLevel = 3;
	size_t orderOfFiniteElement = 2;
	const double dqScaling = 256;

	// chose U (the velocity of the top wall) so that Ma = 0.05
	//const double U = 5. / 100. * sqrt(3) * dqScaling;
	const double U = 5. / 100. * sqrt(3);
	// chose viscosity so that Re = 2000
	//const double Re = 2000;
	const double Re = 2000;
	const double viscosity = U / Re;
	const double startTime = 1;
	const double timeStepSize = 5e-6;
	//const double timeStepSize = 1e-4;

	cout << "Mach number: " << U / (sqrt(3) * dqScaling) << endl;
	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("../results/step-2-stencilbig");
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(1000);
	configuration->setOutputSolutionInterval(500);
	configuration->setNumberOfTimeSteps(100000000);
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
