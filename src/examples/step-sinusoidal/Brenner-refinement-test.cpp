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

	cout << "Starting analysis of sinusoidal shear flow ...." << endl;
	cout
			<< "Usage: compare-Brenner <configuration id> <Ma-Number> <Gamma> <Refinement level>"
			<< endl;

	double cFL = 10;
	bool automatic_decrease = true;
	if (argc != 6) {
		assert(argc == 7);
		cFL = atof(argv[6]);
		automatic_decrease = false;
	}

	const size_t cfg = atoi(argv[1]);
	const double Ma = atof(argv[2]);
	const double gamma = atof(argv[3]);
	const size_t refinementLevel = atoi(argv[4]);
	const size_t orderOfFiniteElement = atoi(argv[5]);
	cout << "Ma = " << Ma << ", gamma = " << gamma << endl;
	const double Re = 1;

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
	fName << getenv("OUTPUT_DIR") << "/result.txt";
	cout << fName.str().c_str();
	std::ofstream resultFile(fName.str().c_str());
	resultFile
			<< "#<configuration id> <Ma-Number> <Gamma> <Refinement level>  eps     alpha      Lx       h       a        b     u_mean    sigma  tau   Psi_s"
			<< endl;

	cout << "Starting configuration " << cfg << "..." << endl;

	/// get geometry parameters
	double Lx = configurations[cfg][0];
	double h = configurations[cfg][1];
	double a = configurations[cfg][2];
	double b = configurations[cfg][3];
	assert(b < h);
	double sigma = sqrt(0.5 * b * b);
	double epsilon = (h - b) / sigma;
	double alpha = h / Lx;

	/// create CFD problem
	double cell_aspect_ratio = 0.25;
	const double viscosity = u_a * h / Re;
	shared_ptr<ProblemDescription<2> > sinusFlow =
			make_shared<SinusoidalShear2D>(viscosity, u_a, refinementLevel, Lx,
					h, b, cell_aspect_ratio);
	const double dt = CFDSolverUtilities::calculateTimestep<2>(
			*sinusFlow->getTriangulation(), orderOfFiniteElement, D2Q9(scaling),
			cFL);

	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("OUTPUT_DIR");
	shared_ptr<SolverConfiguration> configuration = make_shared<
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
	double tau = viscosity / (dt * (scaling * scaling) / 3.0);

	// automatic decrease of timestep to get tau = 0.5
	if ((tau < 0.5) and (automatic_decrease)) {
		double new_dt = viscosity / (scaling * scaling / 3.0);
		tau = viscosity / (new_dt * (scaling * scaling) / 3.0);
		configuration->setTimeStepSize(new_dt);
		cout << "Config " << cfg << ": tau too small ("
				<< viscosity / (dt * (scaling * scaling) / 3.0)
				<< "). Automatic decrease of time step size (now CFL = "
				<< cFL * new_dt / dt << " , tau = " << tau << ")." << endl;

	}

	configuration->setTimeIntegrator(OTHER);
	configuration->setDealIntegrator(CRANK_NICOLSON);

	configuration->setConvergenceThreshold(1e-5 * Ma / sqrt(gamma));
	if (gamma < 1 - 1e-5) {
		configuration->setCollisionScheme(BGK_STEADY_STATE);
		configuration->setBGKSteadyStateGamma(gamma);
	}

	// make solver object and run simulation
	CFDSolver<2> solver(configuration, sinusFlow);
	double Psi_s_old = 0.0;
	double err = 100.0;
	double qx = 0.0;
	double Psi_s = 0.0;
	while (err > 1e-4) {
		for (size_t i = 0; i < 100; i++) {
			solver.stream();
			solver.collide();
		}
		const distributed_vector & ux = solver.getVelocity().at(0);
		qx = u_a * h
				- PhysicalProperties<2>::meanVelocityX(ux,
						solver.getAdvectionOperator()) * h;
		Psi_s = 2 * qx / (sigma * u_a) - h / sigma;
		err = fabs(Psi_s - Psi_s_old);
	}
	resultFile << cfg << "  " << Ma << "  " << refinementLevel << "  " << gamma << "  "  << epsilon << "   " << alpha
			<< "   " << Lx << "   " << h << "   " << a << "   " << b << "   "
			<< qx / h << "   " << sigma << "   " << tau << "  " << Psi_s
			<< endl;
	cout << "Flow factor " << Psi_s << endl;
	cout << "Analysis finished." << endl;

	return 0;
}
