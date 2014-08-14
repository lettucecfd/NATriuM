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
#include "utilities/CFDSolverUtilities.h"

#include "CouetteFlow2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-2..." << endl;


	// set Reynolds and Mach number
	const double Re = 2000;
	const double Ma = 0.1;
	const double tau = 0.9;
	const double dt = 1e-2;

	// set spatial discretization
	size_t refinementLevel = 4;
	size_t orderOfFiniteElement = 2;

	// set Problem so that the right Re and Ma are achieved
	double L = 1.0;
	shared_ptr<CouetteFlow2D> helper = make_shared<CouetteFlow2D>(
			1, 1, refinementLevel, L, 0);
	double dx_min = CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(*(helper->getTriangulation()),orderOfFiniteElement);
	double CFL = sqrt(3) * Ma / Re * L / (dx_min * tau);
	cout << "CFL: " << CFL << endl;
	double dqScaling = CFL * dx_min / dt;
	double U = dqScaling * Ma / sqrt(3);

	//const double U = 1/sqrt(3)*Ma;
	//const double dqScaling = 1;
	const double viscosity = U * L / Re; // (because L = 1)

	// in order to start from a continuous solution, do not start at t=0
	const double startTime = 0.0;

	// set small time step size
	//const double timeStepSize = 0.0001;

	shared_ptr<CouetteFlow2D> couetteFlow = make_shared<CouetteFlow2D>(
			viscosity, U, refinementLevel, 1.0, startTime);
	shared_ptr<Benchmark<2> > couetteProblem = couetteFlow;

	//cout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;
	cout << "CFL: " << dt*dqScaling/CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(*(couetteProblem->getTriangulation()),orderOfFiniteElement);

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
	configuration->setNumberOfTimeSteps(40./dt);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(dt);
	configuration->setCommandLineVerbosity(7);

	BenchmarkCFDSolver<2> solver(configuration, couetteProblem);

	solver.run();

	cout << "NATriuM step-2 terminated." << endl;

	return 0;
}
