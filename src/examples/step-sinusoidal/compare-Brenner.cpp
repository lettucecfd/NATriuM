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

#include "deal.II/numerics/data_out.h"

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"
#include "utilities/CFDSolverUtilities.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "SinusoidalShear2D.h"

using namespace natrium;

// Main function
int main(int argc, char* argv[]) {

	cout << "Starting analysis of sinusoidal shear flow ..." << endl;

	const double Ma = 0.05;
	const double Re = 1.0;
	const double cFL = 10.0;
	const double refinementLevel = 5;
	const double orderOfFiniteElement = 5;
	const double tmax = 5.0;

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
	const size_t n_cfg = 5;
	double configurations[n_cfg][4] = { { 1, 0.3, 0, 0.29 }, { 1, 0.3, 0, 0.1 },
			{ 5, 0.3, 0, 0.1 }, { 1, 0.3, 0, 0.05 }, { 5, 0.3, 0, 0.05 } };

	std::stringstream fName;
	fName << getenv("NATRIUM_HOME") << "/sinus-shear-Brenner1/result.txt";
	std::ofstream resultFile(fName.str().c_str());
	resultFile
			<< "#no     eps     alpha      Lx       h       a        b      Psi_s"
			<< endl;

	/// for all configurations
	for (size_t i; i < n_cfg; i++) {
		cout << "Starting configuration " << i << "..." << endl;

		/// get geometry parameters
		double Lx = configurations[i][0];
		double h = configurations[i][1];
		double a = configurations[i][2];
		double b = configurations[i][3];
		assert(b < h);
		double epsilon = (h - b) / (0.5 * b * b);
		double alpha = h / Lx;

		/// create CFD problem
		const double viscosity = u_a * h / Re;
		shared_ptr<ProblemDescription<2> > sinusFlow = make_shared<
				SinusoidalShear2D>(viscosity, u_a, refinementLevel, Lx, h, b);
		const double dt = CFDSolverUtilities::calculateTimestep<2>(
				*sinusFlow->getTriangulation(), orderOfFiniteElement,
				D2Q9IncompressibleModel(scaling), cFL);

		/// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/sinus-shear-Brenner1/cfg-" << i;
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(10);
		configuration->setOutputCheckpointInterval(1000);
		configuration->setOutputSolutionInterval(100);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(ALL);
		configuration->setTimeStepSize(dt);
		configuration->setDealIntegrator(SDIRK_TWO_STAGES);
		configuration->setTimeIntegrator(OTHER);

		configuration->setInitializationScheme(ITERATIVE);
		configuration->setIterativeInitializationNumberOfIterations(10);
		configuration->setIterativeInitializationResidual(1e-15);

		configuration->setNumberOfTimeSteps(tmax / dt);

		// make solver object and run simulation
		CFDSolver<2> solver(configuration, sinusFlow);
		solver.run();

		// calculate and put out flow factors
		double Psi_s = 0.0;
		const distributed_vector & ux = solver.getVelocity().at(0);
		for (size_t j = 0; j < ux.size(); j++) {
			Psi_s += ux(j);
		}
		// compute average (divide by number of grid points)
		Psi_s /= ux.size();
		// relate to shear velocity
		Psi_s /= u_a;
		resultFile << i << "  " << epsilon << "   " << alpha << "   " << Lx
				<< "   " << h << "   " << a << "   " << b << "   " << Psi_s
				<< endl;
	}
	cout << "Analysis finished." << endl;

	return 0;
}
