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

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-1 ..." << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	const double viscosity = 1;
	// C-E-approach: constant stencil scaling
	// specify Mach number
	const double Ma = 0.05;
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = 4;

	// chose scaling so that the right Ma-number is achieved
	double scaling = sqrt(3) * 1 / Ma;

	const double refinementLevel = 4;
		cout << "refinement Level = " << refinementLevel << endl;
//		for (size_t orderOfFiniteElement = 2; orderOfFiniteElement < 7;
//				orderOfFiniteElement++) {
//			cout << "FE order = " << orderOfFiniteElement << endl;
		// the scaling has to be orders of magnitude greater than the boundary velocity

		// calculate distance between quadrature nodes
		//double dx = 2 * 3.1415926
		//		/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));
		// chose dt so that courant (advection) = 1 for the diagonal directions
		//double dt = dx / (scaling * sqrt(2));
		double dt = 0.001;

		cout << "dt = " << dt << " ...";

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/step-1-iterative";
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(10);
		configuration->setOutputCheckpointInterval(1000);
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

		configuration->setNumberOfTimeSteps(1.0 / dt);

		// make problem and solver objects
		shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel);
		shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
		timestart = clock();
		BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
		time1 = clock() - timestart;

		try {
			solver.run();
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;
			cout << " OK ... Init: " << time1 << " sec; Run: " << time2
					<< " sec." << endl;
		} catch (std::exception& e) {
			cout << " Error" << endl;
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


	cout << "step-1 terminated." << endl;

	return 0;
}
