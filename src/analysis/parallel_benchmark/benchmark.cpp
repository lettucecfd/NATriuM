/**
 * @file benchmark.cpp
 * @short Couette Flow in 2D with time measurement
 * @date 16.11.2015
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

	assert(argc >= 3);

	// set spatial discretization
	size_t refinementLevel = std::atoi(argv[1]);
	size_t orderOfFiniteElement = std::atoi(argv[2]);
	size_t nof_iterations = 200;
	if (argc >= 4){
		nof_iterations = std::atoi(argv[3]);
	}
	if (argc == 5){
		// TODO alternative integrators
	}

	pout << "Performance analysis with N=" << refinementLevel << " and p="
			<< orderOfFiniteElement << endl;

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

	double lups = solver.getNumberOfDoFs() * nof_iterations / (time3 / 1000.0);
	pout
			<< "----------------------------------------------------------------------------------"
			<< endl;
	pout << "all times in ms:" << endl;
	pout
			<< "1)n_mpi_proc  2)N  3)p   4)n_dofs  5)t_build_problem  6)t_build_solver  7)t_per_iteration   8)t_total  9)LUPS  10)LUPS/node"
			<< endl;
	pout << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << " "
			<< refinementLevel << " " << orderOfFiniteElement << " "
			<< solver.getNumberOfDoFs() << time1 << " " << time2 << " "
			<< time3 / nof_iterations << " " << clock() - timestart << lups
			<< " "
			<< lups / dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
			<< endl;
	pout << "done." << endl;

	return 0;
}
