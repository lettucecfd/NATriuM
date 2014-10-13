/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "Cylinder2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-9..." << endl;

	// set Reynolds and Mach number
	const double Re = 10;
	const double Ma = 0.1;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1/sqrt(3)*Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)


	// set small time step size
	const double timeStepSize = 0.01;
	const size_t orderOfFiniteElement = 5;

	cout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;

	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-9";
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setNumberOfTimeSteps(5000./timeStepSize);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(timeStepSize);
	configuration->setCommandLineVerbosity(7);
	//configuration->setDistributionInitType(Iterative);


	shared_ptr<Cylinder2D> cylinder = make_shared<Cylinder2D>(
				viscosity, U);

	shared_ptr<ProblemDescription<2> > couetteProblem = cylinder;
	CFDSolver<2> solver(configuration, couetteProblem);

	solver.run();



	cout << "NATriuM step-9 terminated." << endl;

	return 0;
}
