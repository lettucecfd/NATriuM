/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>
#include <mpi.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/stencils/D2Q9.h"

#include "natrium/utilities/BasicNames.h"

#include "Cylinder2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	cout << "Starting NATriuM step-9..." << endl;
#ifdef WITH_TRILINOS
	static	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
#endif

	// set Reynolds and Mach number
	const double Re = 10;
	const double Ma = 0.1;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1/sqrt(3)*Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)

	// load grid
	shared_ptr<Cylinder2D> cylinder = make_shared<Cylinder2D>(
				viscosity, U);
	D2Q9 stencil(dqScaling);
	// set FE order and time step size
	const size_t orderOfFiniteElement = 2;
	const double cfl=1.5;

	const double timeStepSize = CFDSolverUtilities::calculateTimestep<2>(*cylinder->getMesh(),
			orderOfFiniteElement, stencil, cfl);


	cout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;

	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-9";
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setNumberOfTimeSteps(5000./timeStepSize);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(timeStepSize);
	configuration->setCommandLineVerbosity(7);
	//configuration->setDistributionInitType(Iterative);


	shared_ptr<ProblemDescription<2> > couetteProblem = cylinder;
	CFDSolver<2> solver(configuration, couetteProblem);

	solver.run();



	cout << "NATriuM step-9 terminated." << endl;

	return 0;
}
