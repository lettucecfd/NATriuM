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

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/stencils/D2Q9.h"

#include "LubricationSine.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting analysis of sinusoidal shear flow ..." << endl;

	const double Ma = 0.05;
	const double cFL = 10;

	const double LPhysToSim = 1000; // [1/m]
	const double TPhysToSim = 1e6; // [1/s]
	const double RhoPhysToSim = 1./850; // [m**3/kg]

	// Problem definition
	const double rho = 850 * RhoPhysToSim;
	assert ( fabs(rho-1.0) < 1e-5);
	const double radius = 0.005 * LPhysToSim;
	const double delta_radius = 0.004 * radius; //
	double epsilon = 0.5; // clearance (relative to delta_radius)
	double roughnessHeight = 0;//1e-6 * LPhysToSim;
	size_t roughnessLengthRatio = radius / 2e-6 / LPhysToSim ;
	/*if (argc > 2) {
		epsilon = atof(argv[1]);
	}*/
	const double bottomVelocity = 50.0 * LPhysToSim / TPhysToSim;

	const double amplitude = epsilon * delta_radius;
	const double L = 8 * atan(1) * radius;
	const double eta = 0.032; // dynamic viscosity in Pas
	const double viscosity = eta / (rho / RhoPhysToSim) * LPhysToSim * LPhysToSim / TPhysToSim;
	const double scaling = sqrt(3) * bottomVelocity / Ma;
	const double Re = bottomVelocity * delta_radius / viscosity;
	const double Re_phys = 50.0 * 0.004 * 0.005 / (0.032/850);
	pout << "Physical and dimensionless Reynolds number should agree: " << Re_phys << " " << Re << endl ;
	pout << "epsilon = "  << epsilon << endl;
	// calculation of pressure from simulated density: p - p0 = ( rho - rho0) cs**2 /RhoPhysToSim * TPhysToSim**2 / LPhysToSim**2
	// Solver Definition
	const double refinementLevel = 2;
	const double orderOfFiniteElement = 4;
	const double cellAspectRatio = 5;

	boost::shared_ptr<ProblemDescription<2> > sinusFlow =
			boost::make_shared<LubricationSine>(viscosity, bottomVelocity,
					refinementLevel, L, delta_radius, amplitude, cellAspectRatio, roughnessHeight, roughnessLengthRatio);
	const double dt = CFDSolverUtilities::calculateTimestep<2>(
			*sinusFlow->getMesh(), orderOfFiniteElement,
			D2Q9(scaling), cFL);

	// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/lubrication-sine-0.5-smoothx";
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(1000);
	configuration->setOutputCheckpointInterval(1e9);
	configuration->setOutputSolutionInterval(10000);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setTimeIntegrator(THETA_METHOD);
	configuration->setThetaMethodTheta(0.5);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setTimeStepSize(dt);

	//configuration->setInitializationScheme(ITERATIVE);
	//configuration->setIterativeInitializationNumberOfIterations(100);
	//configuration->setIterativeInitializationResidual(1e-15);

	if (dt > 0.1) {
		pout << "Timestep too big." << endl;
	}

	configuration->setNumberOfTimeSteps(80000000);

	// make solver object
	CFDSolver<2> solver(configuration, sinusFlow);

	solver.run();
	pout << "Analysis finished." << endl;

	return 0;
}
