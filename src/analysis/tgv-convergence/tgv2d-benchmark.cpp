/**
 * @file step-1-variousRe.cpp
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

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	//pout << "Starting NATriuM step-1-variousRe ..." << endl;

	CommandLineParser parser(argc, argv);
	parser.setArgument<double>("Re", "Reynolds number", 10);
	parser.setPositionalArgument<int>("ref-level",
			"Refinement level of the computation grid.");
	parser.setFlag("no_init_rho_analytically",
			"Initialize with constant density (1)");
	parser.setFlag("limiter", "Use a limiter in the advection step");
    parser.setFlag("fneq", "Initialize fneq by gradients");
	parser.setArgument<double>("refine-tol",
			"Refinement tolerance for step size control (only for adaptive time integrators)",
			1e-7);
	parser.setArgument<double>("coarsen-tol",
			"Coarsening tolerance for step size control (only for adaptive time integrators)",
			1e-8);
	parser.setArgument<double>("Ma", "Mach number", 0.05);
    parser.setArgument<double>("horizontal", "U for the horizontal velocity", 0);
	try {
		parser.importOptions();
	} catch (HelpMessageStop&) {
		return 0;
	}

	const double N = parser.getArgument<int>("ref-level");
	const double Ma = parser.getArgument<double>("Ma");
	const double Re = parser.getArgument<double>("Re");
	const double init_rho_analytically = not parser.hasArgument(
			"no_init_rho_analytically");
	int limiter = parser.hasArgument("limiter");
	double refine_tol = parser.getArgument<double>("refine-tol");
	double coarsen_tol = parser.getArgument<double>("coarsen-tol");
	double u_horizontal = parser.getArgument<double>("horizontal");

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////
	const double L = 2 * M_PI;
	const double U = 1;
	const double scaling = sqrt(3) * U / Ma;
	const double viscosity = (L * U) / Re;



    boost::shared_ptr<TaylorGreenVortex2D > tgv_tmp = boost::make_shared<
            TaylorGreenVortex2D>(viscosity, N, U / Ma, init_rho_analytically);
    tgv_tmp->setHorizontalVelocity(u_horizontal);
    boost::shared_ptr<Benchmark<2> > tgv = tgv_tmp;

	//boost::shared_ptr<Benchmark<2> > tgv = boost::make_shared<
	//		TaylorGreenVortex2D>(viscosity, N, U / Ma, init_rho_analytically);


    // setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setRestartAtIteration(0);
	configuration->setUserInteraction(false);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	//configuration->setStencil(Stencil_D2Q25H);
	configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);

	if (limiter) {
		configuration->setVmultLimiter(true);
	}
    if(parser.hasArgument("fneq")){
        configuration->setInitializationScheme(GRADIENTS);
    }


        configuration->setSimulationEndTime(-1.0 / (2.0 * viscosity) * log(0.1));

	parser.applyToSolverConfiguration(*configuration);
	configuration->setEmbeddedDealIntegratorParameters(1.2, 0.8, 0.05,
			configuration->getCFL(), refine_tol, coarsen_tol);

	BenchmarkCFDSolver<2> solver(configuration, tgv);
	double delta_t = solver.getTimeStepSize();

	try {

		double timestart = clock();
		solver.run();
		double runtime = clock() - timestart;
		runtime /= CLOCKS_PER_SEC;
		solver.getErrorStats()->update();
		solver.getSolverStats()->update();
		double kinE_num = solver.getSolverStats()->getKinE();
		//double kinE_ana = M_PI*M_PI*exp(-4*viscosity*solver.getTime());
		double simulated_viscosity = -1.0 / (4 * solver.getTime())
				* log(kinE_num / (M_PI * M_PI));
		double numerical_viscosity = simulated_viscosity - viscosity;
		double u_error = solver.getErrorStats()->getL2VelocityError();
		double rho_error = solver.getErrorStats()->getL2DensityError();
		pout
				<< "N p Ma Re horizontal_u CFL stencil init_rho_analytically  #steps Mean_CFL ||p-p_ana||_inf ||u-u_ana||_2  nu_numerical/nu  runtime"
				<< endl;
		pout << N << " " << configuration->getSedgOrderOfFiniteElement() << " " << Ma << " " << Re << " "
		<< configuration->getTimeIntegrator()
				<< u_horizontal << " "
				<< solver.getConfiguration()->getCFL() << " "
				<< configuration->getStencil() << " "
				<< init_rho_analytically << " " << " " << solver.getIteration()
				<< " "
				<< solver.getTime() / solver.getIteration() / delta_t
						* solver.getConfiguration()->getCFL() << " "
				<< rho_error * (U / Ma) * (U / Ma) << " " << u_error << " "
				<< numerical_viscosity / viscosity << " " << runtime << endl;

	} catch (std::exception& e) {
		pout << " Error" << endl;
	}

	LOG(BASIC) << "NATriuM run complete." << endl;

	//pout << "step-1-variousRe terminated." << endl;

	return 0;
}
