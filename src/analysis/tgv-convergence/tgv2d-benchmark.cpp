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

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	//pout << "Starting NATriuM step-1 ..." << endl;

	const double N = atoi(argv[1]);
	const double p = atoi(argv[2]);
	const double Ma = atof(argv[3]);
	const double Re = atof(argv[4]);
	const double integrator = atoi(argv[5]);
	const double CFL = atof(argv[6]);
	const double collision = atoi(argv[7]);
	const double init_rho_analytically = atoi(argv[8]);
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
	boost::shared_ptr<Benchmark<2> > tgv = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, N, U / Ma, init_rho_analytically);

	double delta_t = CFDSolverUtilities::calculateTimestep<2>(*(tgv->getMesh()),
			p, D2Q9(scaling), CFL);

	// setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setTimeStepSize(delta_t);
	if (collision == 1) {
		configuration->setCollisionScheme(KBC_STANDARD);
	}
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);
	configuration->setEmbeddedDealIntegratorParameters(1.2, 0.8, 0.05 * delta_t,
			delta_t, refine_tol, coarsen_tol);
	// end after Dissipation by one order of magnitude
	// exp(-2vt) = 1/10
	configuration->setSimulationEndTime(-1.0 / (2.0 * viscosity) * log(0.1));
	BenchmarkCFDSolver<2> solver(configuration, tgv);

	try {

		double timestart = clock();
		solver.run();
		double runtime = clock() - timestart;
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
				<< "N p Ma Re integrator CFL collision init_rho_analytically  #steps Mean_CFL ||p-p_ana||_inf ||u-u_ana||_2  nu_numerical/nu  runtime"
				<< endl;
		pout << N << " " << p << " " << Ma << " " << Re << " " << integrator
				<< " " << CFL << " " << collision << " "
				<< init_rho_analytically << " " << " " << solver.getIteration() << " "
				<< solver.getTime() / solver.getIteration()
						/ delta_t * CFL << " "
				<< rho_error * (U / Ma) * (U / Ma) << " " << u_error << " "
				<< numerical_viscosity / viscosity << " " << runtime << endl;

	} catch (std::exception& e) {
		pout << " Error" << endl;
	}

	LOG(BASIC) << "NATriuM run complete." << endl;

	pout << "step-1 terminated." << endl;

	return 0;
}
