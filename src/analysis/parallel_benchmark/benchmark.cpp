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
#include "natrium/stencils/D3Q15.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/CouetteFlow2D.h"
#include "natrium/benchmarks/CouetteFlow3D.h"

#include "natrium/utilities/Info.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Usage: ./benchmark <ref_level> <order_fe> <nof_replicates=1> "
			"<is_3D=false> <nof_iter=200> <integrator_id=1>" << endl;
	assert(argc >= 3);

	// set spatial discretization
	size_t refinementLevel = std::atoi(argv[1]);
	size_t orderOfFiniteElement = std::atoi(argv[2]);
	size_t replicates = 1;
	if (argc >= 4) {
		replicates = std::atoi(argv[3]);
	}
	bool dim_3 = false;
	if (argc >= 5) {
		dim_3 = std::atoi(argv[4]);
		if (!dim_3) {
			pout << "nof_replicates set to 1 in 2D" << endl;
		}
	}
	size_t nof_iterations = 200;
	if (argc >= 6) {
		nof_iterations = std::atoi(argv[5]);
	}
	size_t integrator_id = 1;
	if (argc >= 7) {
		integrator_id = std::atoi(argv[6]);
	}

	// get integrator
	TimeIntegratorName time_integrator;
	DealIntegratorName deal_integrator;
	string integrator_name;
	CFDSolverUtilities::get_integrator_by_id(integrator_id, time_integrator,
			deal_integrator, integrator_name);

	pout << "Performance analysis with N=" << refinementLevel << " and p="
			<< orderOfFiniteElement << endl;
	pout << "Integrator: " << integrator_name << endl;
	pout << "3D?: " << dim_3 << endl;

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
	shared_ptr<ProblemDescription<3> > couetteProblem3D;
	shared_ptr<ProblemDescription<2> > couetteProblem2D;

	// set small time step size
	const double CFL = 0.4;
	double delta_t;
	if (dim_3) {
		couetteProblem3D = make_shared<CouetteFlow3D>(viscosity, U,
				refinementLevel, replicates, startTime, isUnstructured);
		delta_t = CFDSolverUtilities::calculateTimestep<3>(
				*(couetteProblem3D->getMesh()), orderOfFiniteElement,
				D3Q15(dqScaling), CFL);
	} else {
		couetteProblem2D = make_shared<CouetteFlow2D>(viscosity, U,
				refinementLevel, 1.0, startTime, isUnstructured);
		delta_t = CFDSolverUtilities::calculateTimestep<2>(
				*(couetteProblem2D->getMesh()), orderOfFiniteElement,
				D2Q9(dqScaling), CFL);
	}

	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setRestartAtLastCheckpoint(false);
	//configuration->setSwitchOutputOff(true);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setUserInteraction(false);
	configuration->setNumberOfTimeSteps(nof_iterations);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(delta_t);
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);

	size_t n_dofs;
	double lups;
	if (dim_3) {
		configuration->setStencil(Stencil_D3Q15);
		time1 = clock() - timestart;
		CFDSolver<3> solver(configuration, couetteProblem3D);

		// info output
		const vector<dealii::types::global_dof_index>& dofs_per_proc =
				solver.getAdvectionOperator()->getDoFHandler()->n_locally_owned_dofs_per_processor();
		for (size_t i = 0;
				i < dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
				i++) {
			pout << "Process "
					<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
					<< " has " << dofs_per_proc.at(i) << " grid points."
					<< endl;
		}
		cout << "Process "
				<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
				<< " is running on host " << Info::getHostName() << "."
				<< endl;

		time2 = clock() - time1;

		solver.run();
		time3 = clock() - time2;

		lups = solver.getNumberOfDoFs() * nof_iterations / (time3 / 1000.0);
		pout
				<< "----------------------------------------------------------------------------------"
				<< endl;
		pout << "Runtime summary:" << endl;
		solver.printRuntimeSummary();
		n_dofs = solver.getNumberOfDoFs();

	} else {
		time1 = clock() - timestart;
		CFDSolver<2> solver(configuration, couetteProblem2D);

		// info output
		const vector<dealii::types::global_dof_index>& dofs_per_proc =
				solver.getAdvectionOperator()->getDoFHandler()->n_locally_owned_dofs_per_processor();
		for (size_t i = 0;
				i < dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
				i++) {
			pout << "Process "
					<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
					<< " has " << dofs_per_proc.at(i) << " grid points."
					<< endl;
		}
		cout << "Process "
				<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
				<< " is running on host " << Info::getHostName() << "." << endl;

		time2 = clock() - time1;

		solver.run();
		time3 = clock() - time2;

		lups = solver.getNumberOfDoFs() * nof_iterations / (time3 / 1000.0);
		pout
				<< "----------------------------------------------------------------------------------"
				<< endl;
		pout << "Runtime summary:" << endl;
		solver.printRuntimeSummary();
		n_dofs = solver.getNumberOfDoFs();
	}

	// final out
	pout << "all times in ms:" << endl;
	pout
			<< "1)n_mpi_proc  2)N  3)p   4)n_dofs  5)t_build_problem  6)t_build_solver  7)t_per_iteration   8)t_total  9)LUPS  10)LUPS/node"
			<< endl;
	pout << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << " "
			<< refinementLevel << " " << orderOfFiniteElement << " " << n_dofs
			<< " " << time1 << " " << time2 << " " << time3 / nof_iterations
			<< " " << clock() - timestart << " " << lups << " "
			<< lups / dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
			<< endl;
	pout << "done." << endl;

	return 0;
}
