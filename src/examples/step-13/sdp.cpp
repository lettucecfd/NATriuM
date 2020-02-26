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

#include "natrium/solver/BenchmarkCompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/stencils/D2Q9.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/CommandLineParser.h"

#include "SmoothDensityPropagation.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	//pout << "Starting NATriuM step-1 ..." << endl;

	CommandLineParser parser(argc, argv);
	parser.setArgument<double>("Re", "Reynolds number", 10);
	parser.setPositionalArgument<int>("ref-level",
			"Refinement level of the computation grid.");
	parser.setFlag("no_init_rho_analytically",
			"Initialize with constant density (1)");
	parser.setFlag("limiter", "Use a limiter in the advection step");
	parser.setArgument<double>("refine-tol",
			"Refinement tolerance for step size control (only for adaptive time integrators)",
			1e-7);
	parser.setArgument<double>("coarsen-tol",
			"Coarsening tolerance for step size control (only for adaptive time integrators)",
			1e-8);
	parser.setArgument<double>("Ma", "Mach number", 0.05);
    parser.setArgument<double>("horizontal", "U for the horizontal velocity", 0);
    parser.setArgument<double>("visc","viscosity of the fluid",0.000001);
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
	const double L = 1.0;//;2* M_PI;
	const double U = u_horizontal*sqrt(1.4)/sqrt(3.0);
	const double scaling = 1.0;//sqrt(3) * (U)/ Ma;
    double viscosity = parser.getArgument<double>("visc");

    boost::shared_ptr<SmoothDensityPropagation > sdp_tmp = boost::make_shared<
            SmoothDensityPropagation>(viscosity, N, sqrt(1.4)/sqrt(3.0), init_rho_analytically, L, u_horizontal);
    sdp_tmp->setHorizontalVelocity(U);
    boost::shared_ptr<CompressibleBenchmark<2> > sdp = sdp_tmp;

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
	configuration->setStencil(Stencil_D2Q25H);
	//configuration->setOutputSolutionInterval(10);
	configuration->setOutputDirectory("/home/dwilde3m/tmp/sdp");
	//configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);

	if (limiter) {
		configuration->setVmultLimiter(true);
	}

	configuration->setSimulationEndTime((1.0/u_horizontal)/(sqrt(1.4)/sqrt(3.0)));

	parser.applyToSolverConfiguration(*configuration);
	configuration->setEmbeddedDealIntegratorParameters(1.2, 0.8, 0.05,
			configuration->getCFL(), refine_tol, coarsen_tol);

	BenchmarkCompressibleCFDSolver<2> solver(configuration, sdp);
	double delta_t = solver.getTimeStepSize();

	try {

		double timestart = clock();
		solver.run();
		double runtime = clock() - timestart;
		runtime /= CLOCKS_PER_SEC;
		solver.getCompressibleErrorStats()->update();
		solver.getSolverStats()->update();
		double kinE_num = solver.getSolverStats()->getKinE();
		//double kinE_ana = M_PI*M_PI*exp(-4*viscosity*solver.getTime());
		double simulated_viscosity = -1.0 / (4 * solver.getTime())
				* log(kinE_num / (M_PI * M_PI));
		double numerical_viscosity = simulated_viscosity - viscosity;
		double u_error = solver.getCompressibleErrorStats()->getL2VelocityError();
		double rho_error = solver.getCompressibleErrorStats()->getL2DensityError();
		double rho_max = solver.getCompressibleErrorStats()->getMaxDensityError();
		pout
				<< "N p Ma Re horizontal_u CFL stencil init_rho_analytically  #steps Mean_CFL ||p-p_ana||_2 ||p-p_ana||_inf  nu_numerical/nu  runtime"
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
				<< rho_error<< " " << rho_max << " "
				<< numerical_viscosity / viscosity << " " << runtime << endl;

	} catch (std::exception& e) {
		pout << " Error" << endl;
	}

	LOG(BASIC) << "NATriuM run complete." << endl;

	//pout << "step-1 terminated." << endl;

	return 0;
}
