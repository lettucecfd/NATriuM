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

#include "ComplexWall1.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-6..." << endl;

	// set Reynolds and Mach number
	const double Re = 2000;
	const double Ma = 0.1;

	// set spatial discretization
	size_t refinementLevel = 5;
	size_t orderOfFiniteElement = 2;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1/sqrt(3)*Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)

	// set small time step size
	const double timeStepSize = 0.0005;

	cout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;
	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-6";
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setNumberOfTimeSteps(20./timeStepSize);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(timeStepSize);
	configuration->setCommandLineVerbosity(7);
	//configuration->setDistributionInitType(Iterative);

	shared_ptr<ComplexWall1> complex = make_shared<ComplexWall1>(
			viscosity, U, refinementLevel, 1.0);
	shared_ptr<ProblemDescription<2> > prob = complex;
	CFDSolver<2> solver(configuration, prob);

	solver.run();

	cout << "NATriuM step-6 terminated." << endl;

	return 0;
}
