/**
 * @file sinus-shear.cpp
 * @short Sinusoidal Shear flow
 * @date 04.02.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <exception>
#include <ctime>
#include <iomanip>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/solver/PhysicalProperties.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/benchmarks/SinusoidalShear2D.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting analysis of sinusoidal shear flow ...." << endl;
	pout << "Usage: compare-Brenner <Ma-Number> <Gamma> <Refinement level>" << endl;

	double cFL = 10;
	bool automatic_decrease = true;
	if (argc != 4){
		assert (argc == 5);
		cFL = atof(argv[4]);
		automatic_decrease = false;
	}

	const double Ma = atof (argv[1]);
	const double gamma = atof (argv[2]);
	const double refinementLevel = atoi (argv[3]);
	pout << "Ma = " << Ma << ", gamma = " << gamma << endl;
	const double Re = 1;
	const double orderOfFiniteElement = 2 ;

	// parameterization by Brenner:
	// - average height: h
	// - minimal height: delta
	// - amplitude of upper wall: a
	// - amplitude of lower wall: b
	// - Channel length (=wave length): Lx
	// - sigma: standard deviation of wall function, sigma^2 = (a^2 + b^2) / 2
	// - epsilon: relative channel height delta/sigma
	// - alpha: relative channel ratio h/Lx
	// - u_a: velocity of upper wall
	// - u_b: velocity of lower wall
	const double u_a = 0.1;
	const double scaling = sqrt(3) * u_a / Ma;

	// note that we consider only static simulations (u_b = 0, a = 0)
	// with swapped upper and lower wall
	// {Lx, h, a, b}
	const size_t n_cfg = 6;
	// cfg 0 is for calibration (checks flow factor)
	double configurations[n_cfg][4] = { { 1, 0.3, 0, 0 }, { 1, 0.3, 0, 0.1 }, {
			5, 0.3, 0, 0.1 }, { 1, 0.3, 0, 0.29 }, { 1, 0.3, 0, 0.05 }, { 5,
			0.3, 0, 0.05 } };

	std::stringstream fName;
	fName << getenv("NATRIUM_HOME") << "/Brenner/Ma-" << std::fixed << std::setprecision(5) << Ma << "-ref-" << std::fixed << std::setprecision(0) << refinementLevel << "/result.txt";
	pout << fName.str().c_str();
	std::ofstream resultFile(fName.str().c_str());
	resultFile
			<< "#gamma   no    eps     alpha      Lx       h       a        b     u_mean    sigma  tau   Psi_s"
			<< endl;

	/// for all configurations
	for (size_t i = 1; i < n_cfg; i++) {
		// Gamma falsifies the velocity field (flow factor drops)
		//for (double gamma = 0.005; gamma <= 1; gamma += 0.15) {
			//size_t i = 0;
			pout << "Starting configuration " << i << "..." << endl;

			/// get geometry parameters
			double Lx = configurations[i][0];
			double h = configurations[i][1];
			double a = configurations[i][2];
			double b = configurations[i][3];
			assert(b < h);
			double sigma = sqrt(0.5 * b * b);
			double epsilon = (h - b) / sigma;
			double alpha = h / Lx;

			/// create CFD problem
			double cell_aspect_ratio = 0.25;
			const double viscosity = u_a * h / Re;
			boost::shared_ptr<ProblemDescription<2> > sinusFlow = boost::make_shared<
					SinusoidalShear2D>(viscosity, u_a, refinementLevel, Lx, h,
					b, cell_aspect_ratio);
			const double dt = CFDSolverUtilities::calculateTimestep<2>(
					*sinusFlow->getMesh(), orderOfFiniteElement,
					D2Q9(scaling), cFL);

			/// setup configuration
			std::stringstream dirName;
			dirName << getenv("NATRIUM_HOME") << "/Brenner/Ma-" << std::fixed << std::setprecision(5) << Ma << "-ref-" << std::fixed <<  std::setprecision(0) << refinementLevel << "/gamma-"
					<< std::fixed <<  std::setprecision(4) << gamma << "_cfg-" << std::fixed <<  std::setprecision(0) << i;
			boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
					SolverConfiguration>();
			//configuration->setSwitchOutputOff(true);
			configuration->setOutputDirectory(dirName.str());
			configuration->setRestartAtLastCheckpoint(false);
			configuration->setUserInteraction(false);
			configuration->setOutputTableInterval(1);
			configuration->setOutputCheckpointInterval(100000000);
			configuration->setOutputSolutionInterval(1000000);
			configuration->setCommandLineVerbosity(WELCOME);
			configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
			configuration->setStencilScaling(scaling);
			configuration->setCommandLineVerbosity(ALL);
			configuration->setTimeStepSize(dt);
			double tau = viscosity/(dt*(scaling*scaling)/3.0);

			if ((tau < 0.5) and (automatic_decrease)) {
				double new_dt  = viscosity/(scaling*scaling/3.0);
				tau = viscosity/(new_dt*(scaling*scaling)/3.0);
				configuration->setTimeStepSize(new_dt);
				pout << "Config " << i <<": tau too small. Automatic decrease of time step size (now CFL = " << cFL * new_dt / dt << " , tau = " << tau << ")." << endl;

			}



			configuration->setTimeIntegrator(OTHER);
			configuration->setDealIntegrator(CRANK_NICOLSON);

			//configuration->setInitializationScheme(ITERATIVE);
			//configuration->setIterativeInitializationNumberOfIterations(100);
			//configuration->setIterativeInitializationResidual(1e-15);

			configuration->setConvergenceThreshold(1e-5*Ma/sqrt(gamma));
			if (gamma < 1 - 1e-5){
				configuration->setCollisionScheme(BGK_STEADY_STATE);
				configuration->setBGKSteadyStateGamma(gamma);
			}

			// make solver object and run simulation
			CFDSolver<2> solver(configuration, sinusFlow);
			solver.run();

			// calculate and put out flow factors
			const distributed_vector & ux = solver.getVelocity().at(0);
			// flow factor
			double qx = u_a*h - PhysicalProperties<2>::meanVelocityX(ux, solver.getAdvectionOperator())*h;
			double Psi_s = 2*qx / (sigma * u_a) - h / sigma;
			resultFile << gamma << "  " << i << "  " << epsilon << "   "
					<< alpha << "   " << Lx << "   " << h << "   " << a << "   "
					<< b << "   " << qx/h << "   " << sigma << "   " << tau << "  " << Psi_s  << endl;
			pout << "Flow factor " << Psi_s << endl;
		//}
	}
	pout << "Analysis finished." << endl;

	return 0;
}
