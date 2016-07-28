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
#include "natrium/stencils/D3Q19.h"

#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/TaylorGreenVortex3D.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

#include "natrium/utilities/Info.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Usage: ./benchmark <ref_level> <order_fe> "
			"<nof_iter=200>" << endl;
	assert(argc >= 3);

	// set spatial discretization
	size_t refinement_level = std::atoi(argv[1]);
	size_t order_fe = std::atoi(argv[2]);
	size_t nof_iterations = 200;
	if (argc >= 4) {
		nof_iterations = std::atoi(argv[3]);
	}

	pout << "Performance analysis with N=" << refinement_level << " and p="
			<< order_fe << endl;

	// set Reynolds and Mach number
	const double Re = 1;
	const double Ma = 0.05;

	// set Problem so that the right Re and Ma are achieved
	const double U = 1;
	const double L = 2 * M_PI;
	const double viscosity = U * L / Re;
	const double cs = U / Ma;

	// chose scaling so that the right Ma-number is achieved
	const double scaling = sqrt(3) * cs;
	const bool init_rho_analytically = true;

	// time measurement variables
	double time1, time2, time3, timestart;
	timestart = clock();
	boost::shared_ptr<ProblemDescription<3> > tgvProblem3D;

	// set small time step size
	const double CFL = 0.4;
	pout << "Make problem..." << endl;
	tgvProblem3D = boost::make_shared<TaylorGreenVortex3D>(viscosity,
			refinement_level, cs, init_rho_analytically);
	pout << "...done" << endl;

	// configure solver
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setRestartAtLastCheckpoint(false);
	//configuration->setSwitchOutputOff(true);
	//configuration->setCommandLineVerbosity(ALL);
	configuration->setUserInteraction(false);
	configuration->setNumberOfTimeSteps(nof_iterations);
	configuration->setSedgOrderOfFiniteElement(order_fe);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	configuration->setCommandLineVerbosity(DETAILED);

	size_t n_dofs;
	double lups;
	configuration->setStencil(Stencil_D3Q19);
	time1 = clock() - timestart;
	pout << "Make solver..." << endl;
	CFDSolver<3> solver(configuration, tgvProblem3D);
	pout << "...done" << endl;

	// info output
	const vector<dealii::types::global_dof_index>& dofs_per_proc =
			solver.getAdvectionOperator()->getDoFHandler()->n_locally_owned_dofs_per_processor();

	cout << "Process "
			<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
			<< " has "
			<< solver.getAdvectionOperator()->getDoFHandler()->n_locally_owned_dofs()
			<< " grid points." << endl;
	cout << "Process "
			<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
			<< " is running on host " << Info::getHostName() << "." << endl;

	time2 = clock() - time1;

	solver.run();
	time3 = clock() - time2;

	lups = solver.getNumberOfDoFs() * nof_iterations / (time3 /CLOCKS_PER_SEC);
	pout
			<< "----------------------------------------------------------------------------------"
			<< endl;
	pout << "Runtime summary:" << endl;
	solver.printRuntimeSummary();
	n_dofs = solver.getNumberOfDoFs();

	// final out
	pout << "all times in ms:" << endl;
	pout
			<< "1)n_mpi_proc  2)N  3)p   4)#grid points  5)t_build_problem  6)t_build_solver  7)t_per_iteration   8)t_total  9)LUPS  10)LUPS/process"
			<< endl;
	pout << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << " "
			<< refinement_level << " " << order_fe << " " << n_dofs << " "
			<< time1/CLOCKS_PER_SEC*1000 << " " << time2/CLOCKS_PER_SEC*1000 << " " << time3/CLOCKS_PER_SEC*1000 / nof_iterations << " "
			<< (clock() - timestart)/CLOCKS_PER_SEC*1000 << " " << lups << " "
			<< lups / dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
			<< endl;
	pout << "done." << endl;

	return 0;
}
