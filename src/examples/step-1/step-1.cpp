/**
 * @file step-1.cpp
 * @short First tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "TaylorGreenVortex2D.h"

using namespace natrium;


// Main function
int main() {

	cout << "Starting NATriuM step-1..." << endl;

	// set parameters, set up configuration object
	size_t refinementLevel = 5;
	size_t orderOfFiniteElement = 2;
	double viscosity = 1;

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("../results/step-1");
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputTableInterval(100);
	configuration->setOutputCheckpointInterval(100);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(50);
	double tScaling = std::min(0.1, 1. / (2 * configuration->getStencilScaling()));
	double deltaX = 1.
			/ (pow(2, refinementLevel)
					* (configuration->getSedgOrderOfFiniteElement() - 1));
	configuration->setTimeStepSize(tScaling * deltaX);
	configuration->setNumberOfTimeSteps(5000);

	//configuration->setDistributionInitType(Iterative);

	// make problem and solver objects
	shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel);
	shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
	BenchmarkCFDSolver<2> solver(configuration, taylorGreen);


	solver.run();

	cout << "NATriuM step-1 terminated." << endl;

	return 0;
}
