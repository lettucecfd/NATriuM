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

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CommandLineParser.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// parse arguments from command line
	CommandLineParser parse(argc, argv);
	parse.addDocumentationString("step-1-taylorGreenVortex", "This program runs the two-dimensional Taylor-Green vortex on a regular grid.");
	parse.setArgument<double>("L", "size of the computational domain", 2*M_PI);
	parse.setArgument<int>("ref-level", "refinement level of the computational grid", 5);
	try {
		parse.importOptions();
	} catch (HelpMessageStop&){
		return 0;
	}


	pout << "Starting NATriuM step-1 ..." << endl;
	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	//const double viscosity = 1;
	double L = parse.getArgument<double>("L");
	const double Re = 10;
	const double viscosity = 1 * L / Re;
	// C-E-approach: constant stencil scaling
	// specify Mach number
	const double Ma = 0.05;
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = 1;

	// chose scaling so that the right Ma-number is achieved
	double scaling = sqrt(3) * 2 / Ma;

	const double refinementLevel = parse.getArgument<int>("ref-level");
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
		double CFL=sqrt(2)*4;

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/step-1-iterative";
		boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(5);
		configuration->setOutputCheckpointInterval(1000000000);
		configuration->setOutputSolutionInterval(1);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(ALL);
		configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
		configuration->setCFL(CFL);

      //  configuration->setInitializationScheme(ITERATIVE);
        configuration->setIterativeInitializationNumberOfIterations(1000);
        configuration->setIterativeInitializationResidual(1e-7);

		configuration->setSimulationEndTime(10.0);


		parse.applyToSolverConfiguration(*configuration);
		// make problem and solver objects
		boost::shared_ptr<TaylorGreenVortex2D> tgVortex = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel, 1./Ma, true, L);
		//tgVortex->setHorizontalVelocity(1);
		boost::shared_ptr<ProblemDescription<2> > taylorGreen = tgVortex;
		timestart = clock();

		pout << "Make solver" << endl;
		CFDSolver<2> solver(configuration, taylorGreen);
		time1 = clock() - timestart;

		pout << "Check for standard LBM" << endl;
		distributed_block_vector ones;
		distributed_block_vector result;
		ones.reinit(8);
		result.reinit(8);
		for (size_t i = 0; i < 8; i++) {
			ones.block(i).reinit(solver.getAdvectionOperator()->getLocallyOwnedDofs(), MPI_COMM_WORLD);
			result.block(i).reinit(solver.getAdvectionOperator()->getLocallyOwnedDofs(), MPI_COMM_WORLD);
			// reinit does only change the size but not the content
			//for all degrees of freedom on current processor
			dealii::IndexSet::ElementIterator it(solver.getAdvectionOperator()->getLocallyOwnedDofs().begin());
			dealii::IndexSet::ElementIterator end(solver.getAdvectionOperator()->getLocallyOwnedDofs().end());
			for (; it != end; it++) {
				size_t j = *it;
				ones.block(i)(j) = 1;
			}
		}
        //solver.getAdvectionOperator()->getSystemMatrix().print(cout);
        //solver.getAdvectionOperator()->getSystemMatrix().vmult(result, ones);

		result -= ones;
		cout << "error: " << result.norm_sqr();
		//solver.getF().getFStream().print(cout,10,true);


		cout << "Run simulation" << endl;
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
