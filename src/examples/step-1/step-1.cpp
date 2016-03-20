/**
 * @file step-1.cpp
 * @short Simulation of the Taylor-Green decaying vortex in 2D (only periodic walls).
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/Benchmark.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// ========================================================================
	// READ COMMAND LINE PARAMETERS
	// ========================================================================
	pout
			<< "Usage: ./step-1 <refinement_level=3> <p=4> <collision-id=0 (BGK: 0, KBC: 1)> <filter=0 (no: 0, exp: 1, new: 2> <integrator-id=1> <CFL=0.4> <stencil_scaling=1.0>"
			<< endl;

	size_t refinement_level = 2;
	if (argc >= 2) {
		refinement_level = std::atoi(argv[1]);
	}
	pout << "... N:    " << refinement_level << endl;

	size_t p = 2;
	if (argc >= 3) {
		p = std::atoi(argv[2]);
	}
	pout << "... p:    " << p << endl;

	size_t collision_id = 0;
	if (argc >= 4) {
		collision_id = std::atoi(argv[3]);
	}
	pout << "... Coll:  " << collision_id << endl;

	size_t filter_id = 0;
	if (argc >= 5) {
		filter_id = std::atoi(argv[4]);
	}
	pout << "... Filter:  " << filter_id << endl;

	size_t integrator_id = 1;
	if (argc >= 6) {
		integrator_id = std::atoi(argv[5]);
	}
	pout << "... Int:  " << integrator_id << endl;

	double CFL = .1;
	if (argc >= 7) {
		CFL = std::atof(argv[6]);
	}
	pout << "... CFL:  " << CFL << endl;

	double stencil_scaling = 1.0;
	if (argc >= 8) {
		stencil_scaling = std::atof(argv[7]);
	}
	pout << "... stencil_scaling:  " << stencil_scaling << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// get integrator
	TimeIntegratorName time_integrator;
	DealIntegratorName deal_integrator;
	string integrator_name;
	CFDSolverUtilities::get_integrator_by_id(integrator_id, time_integrator,
			deal_integrator, integrator_name);

	// Re = viscosity/(2*pi)
	const double viscosity = 1;
	// C-E-approach: constant stencil scaling
	// specify Mach number
	const double Ma = 0.05;
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = p;

	// chose scaling so that the right Ma-number is achieved
	double scaling = sqrt(3) * 1 / Ma;
	CFL /=scaling;

	const double refinementLevel = refinement_level;
		pout << "refinement Level = " << refinementLevel << endl;
//		for (size_t orderOfFiniteElement = 2; orderOfFiniteElement < 7;
//				orderOfFiniteElement++) {
//			pout << "FE order = " << orderOfFiniteElement << endl;
		// the scaling has to be orders of magnitude greater than the boundary velocity

		// calculate distance between quadrature nodes
		//double dx = 2 * 3.1415926
		//		/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));
		// chose dt so that courant (advection) = 1 for the diagonal directions
		//double dt = dx / (scaling * sqrt(2));


		double dt = 0.001;



		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/step-1-iterative";
		boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(10);
		configuration->setOutputCheckpointInterval(1000);
		configuration->setOutputSolutionInterval(10);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(ALL);
		configuration->setTimeIntegrator(time_integrator);
		configuration->setDealIntegrator(deal_integrator);

		configuration->setCollisionScheme(BGK_STANDARD);
		if (collision_id == 1) {
			configuration->setCollisionScheme(KBC_STANDARD);
		}

		//configuration->setInitializationScheme(ITERATIVE);
		//configuration->setIterativeInitializationNumberOfIterations(1000);
		//configuration->setIterativeInitializationResidual(1e-15);





		// make problem and solver objects
		boost::shared_ptr<TaylorGreenVortex2D> tgVortex = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel, 1./Ma);
		dt = CFDSolverUtilities::calculateTimestep<2>(
					*(tgVortex->getMesh()), p, D2Q9(stencil_scaling), CFL);
		if (dt > 0.1) {
			pout << "Timestep too big." << endl;
		}
		pout << "dt = " << dt << " ...";
		configuration->setTimeStepSize(dt);
		configuration->setNumberOfTimeSteps(10.0 / dt);
		boost::shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
		timestart = clock();
		BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
		time1 = clock() - timestart;

		try {
			solver.run();
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;
			pout << " OK ... Init: " << time1 << " sec; Run: " << time2
					<< " sec." << endl;
		} catch (std::exception& e) {
			pout << " Error" << endl;
		}

/*
		size_t N = solver.getConfiguration()->getNumberOfTimeSteps();
		for (size_t i = 0; i < N; i++) {
			solver.output(i);
			solver.stream();
			solver.collide();
			if (i % 10 == 0) {
						std::stringstream str;
						str << dirName.str().c_str() << "/f_"
								<< i << ".vtu";
						std::string filename = str.str();
						std::ofstream vtu_output(filename.c_str());
						dealii::DataOut<2> data_out;
						data_out.attach_dof_handler(*solver.getAdvectionOperator()->getDoFHandler());
						data_out.add_data_vector(solver.getF().at(1), "f1");
						data_out.add_data_vector(solver.getF().at(2), "f2");
						data_out.add_data_vector(solver.getF().at(3), "f3");
						data_out.add_data_vector(solver.getF().at(4), "f4");
						data_out.add_data_vector(solver.getF().at(5), "f5");
						data_out.add_data_vector(solver.getF().at(6), "f6");
						data_out.add_data_vector(solver.getF().at(7), "f7");
						data_out.add_data_vector(solver.getF().at(8), "f8");

						data_out.build_patches(20);
						data_out.write_vtu(vtu_output);
					}
		}
		solver.output(N);*/
		LOG(BASIC) << "NATriuM run complete." << endl;


	pout << "step-1 terminated." << endl;

	return 0;
}
