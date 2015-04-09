/**
 * @file sinus-shear.cpp
 * @short Sinusoidal Shear flow
 * @date 04.02.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
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

#include "../../examples/step-sinusoidal/SinusoidalShear2D.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	cout << "Starting analysis of sinusoidal shear flow ..." << endl;

	const double Ma = 0.1;
	const double cFL = 1.0;
	const double alpha = 0.2;
	const double epsilon = 0.2;
	const double Re = 1.0 * 2 / epsilon / alpha;  // With reduced Re 1.0 (Hu, 1997): Re = Re* * 2 / epsilon / alpha;
	const double bottomVelocity = 1.0;
	const double averageHeight = 1.0;
	const double amplitude = epsilon * averageHeight;
	const double L = amplitude/alpha;
	const double viscosity = 2 * bottomVelocity * L / Re;
	const double scaling = bottomVelocity / Ma;
	const double refinementLevel = 4;
	const double orderOfFiniteElement = 5;

	shared_ptr<ProblemDescription<2> > sinusFlow =
			make_shared<SinusoidalShear2D>(viscosity, bottomVelocity,
					refinementLevel, L, averageHeight, amplitude);
	const double dt = CFDSolverUtilities::calculateTimestep<2>(
			*sinusFlow->getTriangulation(), orderOfFiniteElement,
			D2Q9(scaling), cFL);

	// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/sinus-shear";
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

	configuration->setNumberOfTimeSteps(10.0 / dt);

	// make solver object
	CFDSolver<2> solver(configuration, sinusFlow);

	solver.run();
	cout << "Analysis finished." << endl;

	return 0;
}
