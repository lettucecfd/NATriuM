/*
 * couette-convergence.cpp
 *
 *  Created on: 18.02.2018
 *      Author: akraemer
 */




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
#include "natrium/stencils/D2Q9.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/CommandLineParser.h"

#include "natrium/benchmarks/CouetteFlow2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	//pout << "Starting NATriuM step-1 ..." << endl;

	CommandLineParser parser(argc, argv);
	parser.setArgument<double>("Re", "Reynolds number", 2000);
	parser.setPositionalArgument<int>("ref-level",
			"Refinement level of the computation grid.");
	parser.setArgument<double>("Ma", "Mach number", 0.05);
	try {
		parser.importOptions();
	} catch (HelpMessageStop&) {
		return 0;
	}

	const double N = parser.getArgument<int>("ref-level");
	const double Ma = parser.getArgument<double>("Ma");
	const double Re = parser.getArgument<double>("Re");
	const double start_time = 0.0;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////
	const double L = 1;
	const double U = 1;
	const double scaling = sqrt(3) * U / Ma;
	const double viscosity = (L * U) / Re;

	boost::shared_ptr<Benchmark<2> > couette = boost::make_shared<
			CouetteFlow2D>(viscosity, U, N, L, start_time);

	// setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setRestartAtIteration(0);
	configuration->setUserInteraction(false);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);


	configuration->setSimulationEndTime(40.0);

	parser.applyToSolverConfiguration(*configuration);

	BenchmarkCFDSolver<2> solver(configuration, couette);
	double delta_t = solver.getTimeStepSize();

	try {

		double timestart = clock();
		solver.run();
		double runtime = clock() - timestart;
		runtime /= CLOCKS_PER_SEC;
		solver.getErrorStats()->update();
		solver.getSolverStats()->update();
		// double kinE_num = solver.getSolverStats()->getKinE();
		//double kinE_ana = M_PI*M_PI*exp(-4*viscosity*solver.getTime());
		double u_error = solver.getErrorStats()->getL2VelocityError();
		double rho_error = solver.getErrorStats()->getL2DensityError();
		pout
				<< "N p Ma Re integrator CFL collision  #steps Mean_CFL ||p-p_ana||_inf ||u-u_ana||_2  runtime"
				<< endl;
		pout << N << " " << configuration->getSedgOrderOfFiniteElement() << " " << Ma << " " << Re << " "
				<< "td" << configuration->getTimeIntegrator()
				<< configuration->getDealIntegrator() << " "
				<< solver.getConfiguration()->getCFL() << " "
				<< configuration->getCollisionScheme() << " "
				<< " " << solver.getIteration()
				<< " "
				<< solver.getTime() / solver.getIteration() / delta_t
						* solver.getConfiguration()->getCFL() << " "
				<< rho_error * (U / Ma) * (U / Ma) << " " << u_error << " "
			    << runtime << endl;

	} catch (std::exception& e) {
		pout << " Error" << endl;
	}

	LOG(BASIC) << "NATriuM run complete." << endl;

	//pout << "step-1 terminated." << endl;

	return 0;
}
