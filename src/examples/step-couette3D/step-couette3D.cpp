/**
 * @file step-couette3D.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/stencils/D3Q19.h"

#include "natrium/benchmarks/CouetteFlow3D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-couette3D..." << endl;

	// set Reynolds and Mach number
	const double Re = 550;
	const double Ma = 0.1;

	// set spatial discretization
	size_t refinementLevel = 3;
	size_t orderOfFiniteElement = 4;
	bool isUnstructured = false;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1/sqrt(3)*Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)

	// in order to start from a continuous solution, do not start at t=0
	const double startTime = 0.0;

	boost::shared_ptr<Benchmark<3> > couetteProblem = boost::make_shared<CouetteFlow3D>(
			viscosity, U, refinementLevel, 1.0, startTime, isUnstructured);

	// set small time step size
	const double CFL = 0.4;

	pout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;
	// configure solver
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-couette3D";
	configuration->setStencil(Stencil_D3Q19);
	configuration->setOutputDirectory(dirname.str());
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setNumberOfTimeSteps(1000);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setCFL(CFL);
	configuration->setCommandLineVerbosity(7);
	//configuration->setDistributionInitType(Iterative);

	BenchmarkCFDSolver<3> solver(configuration, couetteProblem);

	solver.run();

	pout << "NATriuM step-2 terminated." << endl;

	return 0;
}
