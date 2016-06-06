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
#include "deal.II/grid/grid_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/stencils/D2Q9.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "CompareToBotella.h"
#include "../../examples/step-0/LidDrivenCavity2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	//pout << "Starting NATriuM  ..." << endl;

	const double N = atoi(argv[1]);
	const double p = atoi(argv[2]);
	const double Ma = atof(argv[3]);
	const size_t Re = atoi(argv[4]);
	const double integrator = atoi(argv[5]);
	const double CFL = atof(argv[6]);
	const double collision = atoi(argv[7]);
	const bool semi_lagrange = atoi(argv[8]);
	double refine_tol = 1e-7;
	double coarsen_tol = 1e-8;
	if (argc > 9) {
		refine_tol = atof(argv[9]);
	}
	if (argc > 10) {
		refine_tol = atof(argv[10]);
	}

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////
	const double L = 2 * M_PI;
	const double U = 1;
	const double scaling = sqrt(3) * U / Ma;
	const double viscosity = (L * U) / Re;
	TimeIntegratorName time_integrator;
	DealIntegratorName deal_integrator;
	string integrator_name;
	CFDSolverUtilities::get_integrator_by_id(integrator, time_integrator,
			deal_integrator, integrator_name);
	pout << "... that is the " << integrator_name << endl;
	pout << "-------------------------------------" << endl;
	boost::shared_ptr<ProblemDescription<2> > tgv = boost::make_shared<
			LidDrivenCavity2D>(U, viscosity, N);
	double delta_t = CFDSolverUtilities::calculateTimestep<2>(*(tgv->getMesh()),
			p, D2Q9(scaling), CFL);


	// setup configuration
	std::stringstream outdir;
	outdir << getenv("NATRIUM_HOME") << "/cavity-benchmark/Re" << Re << "_N"
			<< N << "_p" << p << "_Ma" << Ma << "_int" << integrator << "_CFL"
			<< CFL << "_coll" << collision << "_sl" << semi_lagrange;
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory(outdir.str());
	configuration->setOutputSolutionInterval(100000000);
	configuration->setOutputCheckpointInterval(100000);
	configuration->setOutputTableInterval(1000);
	//configuration->setSwitchOutputOff(true);
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setCommandLineVerbosity(BASIC);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setCFL(CFL);
	if (collision == 1) {
		configuration->setCollisionScheme(KBC_STANDARD);
	}
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);
	configuration->setEmbeddedDealIntegratorParameters(1.2, 0.8, 0.05,
			CFL, refine_tol, coarsen_tol);
	// end after Dissipation by one order of magnitude
	// exp(-2vt) = 1/10
	//configuration->setSimulationEndTime(-1.0 / (2.0 * viscosity) * log(0.1));
	configuration->setSimulationEndTime(50);
	configuration->setConvergenceThreshold(1e-10);
	if (semi_lagrange){
		configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	}
	CFDSolver<2> solver(configuration, tgv);
	std::stringstream outfile;
	outfile << getenv("NATRIUM_HOME") << "/cavity-benchmark/Re" << Re << "_N"
			<< N << "_p" << p << "_Ma" << Ma << "_int" << integrator << "_CFL"
			<< CFL << "_coll" << collision << "_sl" << semi_lagrange << ".dat";
	boost::shared_ptr<CompareToBotella> res = boost::make_shared<CompareToBotella>(
			solver, Re, outfile.str());
	solver.appendDataProcessor(res);
	// put out grid
	std::stringstream gridfile;
	gridfile << getenv("NATRIUM_HOME") << "/cavity-benchmark/grid_N"
			<< N << "_p" << p << ".eps";
	std::ofstream gout(gridfile.str());
	dealii::GridOut grid_out;
	grid_out.write_eps(*tgv->getMesh(), gout);
	gout.close();
	// run simulation
	try {

		double timestart = clock();
		solver.run();
		double runtime = clock() - timestart;
		solver.getSolverStats()->update();
		res->apply();
		res->printFinalVelocities();
		double u_error = res->getUError();
		double v_error = res->getVError();
		pout
				<< "N p Ma Re integrator CFL collision  #steps Mean_CFL u_error v_error  runtime"
				<< endl;
		pout << N << " " << p << " " << Ma << " " << Re << " " << integrator
				<< " " << CFL << " " << collision << " "
				<< solver.getIteration() << " "
				<< solver.getTime() / solver.getIteration() / delta_t * CFL
				<< " " << u_error << " " << v_error << " " << " " << runtime
				<< endl;

	} catch (std::exception& e) {
		pout << " Error" << endl;
	}

	LOG(BASIC) << "NATriuM run complete." << endl;


	return 0;
}
