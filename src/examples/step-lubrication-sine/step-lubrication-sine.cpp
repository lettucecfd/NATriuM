/**
 * @file step-lubrication-sine.cpp
 * @short Sinusoidal Shear flow
 * @date 10.02.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <exception>
#include <ctime>

#include "deal.II/numerics/data_out.h"

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"
#include "utilities/CFDSolverUtilities.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "LubricationSine.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	cout << "Starting analysis of sinusoidal shear flow ..." << endl;

	const double Ma = 0.1;
	const double cFL = 0.4;
	const double factorL = 1000;
	const double factorT = 1e6;
	const double factorRho = 1./850;

	// Problem definition
	const double radius = 0.005 * factorL;
	const double delta_radius = 0.004 * radius; //
	double epsilon = 0.5; // clearance (relative to delta_radius)
	/*if (argc > 2) {
		epsilon = atof(argv[1]);
	}*/
	const double bottomVelocity = 50.0 * factorL / factorT;

	const double amplitude = epsilon * delta_radius;
	const double L = 8 * atan(1) * radius;
	const double eta = 0.032; // dynamic viscosity in Pas
	const double density = 1;
	const double viscosity = eta * factorRho / density;
	const double scaling = bottomVelocity / Ma;
	const double Re = bottomVelocity * L / viscosity;
	cout << "epsilon = "  << epsilon << endl;

	// Solver Definition
	const double refinementLevel = 1;
	const double orderOfFiniteElement = 1;
	const double cellAspectRatio = 10;

	shared_ptr<ProblemDescription<2> > sinusFlow =
			make_shared<LubricationSine>(viscosity, bottomVelocity,
					refinementLevel, L, delta_radius, amplitude, cellAspectRatio);
	const double dt = CFDSolverUtilities::calculateTimestep<2>(
			*sinusFlow->getTriangulation(), orderOfFiniteElement,
			D2Q9IncompressibleModel(scaling), cFL);

	// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/lubrication-sine-0.5";
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(10);
	configuration->setOutputCheckpointInterval(1000);
	configuration->setOutputSolutionInterval(100);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setTimeStepSize(dt);

	configuration->setInitializationScheme(ITERATIVE);
	configuration->setIterativeInitializationNumberOfIterations(1000);
	configuration->setIterativeInitializationResidual(1e-15);

	if (dt > 0.1) {
		cout << "Timestep too big." << endl;
	}

	configuration->setNumberOfTimeSteps(800000);

	// make solver object
	CFDSolver<2> solver(configuration, sinusFlow);

	solver.run();
	cout << "Analysis finished." << endl;

	return 0;
}
