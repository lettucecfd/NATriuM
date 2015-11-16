/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"
#include "deal.II/base/utilities.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/CouetteFlow2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	assert(argc == 3);

	// set spatial discretization
	size_t refinementLevel = std::atoi(argv[1]);
	size_t orderOfFiniteElement = std::atoi(argv[2]);

	pout << "Performance analysis with N=" << refinementLevel << " and p="
			<< orderOfFiniteElement << endl;
	pout << "on " << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
			<< " mpi processes:" << endl;

	bool isUnstructured = false;

	// set Reynolds and Mach number
	const double Re = 2000;
	const double Ma = 0.05;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1 / sqrt(3) * Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)

	// in order to start from a continuous solution, do not start at t=0
	const double startTime = 0.0;
	const double nof_iterations = 50;

	// time measurement variables
	double time1, time2, time3, timestart;
	timestart = clock();
	shared_ptr<ProblemDescription<2> > couetteProblem = make_shared<
			CouetteFlow2D>(viscosity, U, refinementLevel, 1.0, startTime,
			isUnstructured);

	// set small time step size
	const double CFL = 0.4;
	const double timeStepSize = CFDSolverUtilities::calculateTimestep<2>(
			*(couetteProblem->getMesh()), orderOfFiniteElement, D2Q9(dqScaling),
			CFL);

	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-2";
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setSwitchOutputOff(true);
	configuration->setUserInteraction(false);
	configuration->setNumberOfTimeSteps(nof_iterations);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(timeStepSize);

	time1 = clock() - timestart;
	CFDSolver<2> solver(configuration, couetteProblem);
	time2 = clock() - time1;

	solver.run();
	time3 = clock() - time2;

	double mlups = 1e-6 * solver.getNumberOfDoFs() / time3;
	pout << "----------------------------------------------------------------------------------"
			<< endl;
	pout << "all times in ms:" << endl;
	pout << "n_mpi_proc         t_build_solver         t_per_iteration      MLUPS       t_total" << endl;
	pout << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << " "
			<< time1 << " " << time2 << " " << time3 / nof_iterations << " " << mlups << " "
			<< clock() - timestart << endl;
	pout << "done." << endl;

	return 0;
}