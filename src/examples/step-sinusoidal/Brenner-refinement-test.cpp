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
	pout
			<< "Usage: compare-Brenner <configuration id> <Ma-Number> <Gamma> <Refinement level> <order FE> [<CFL>] [<Integrator ID>]"
			<< endl <<
			"0: FORWARD_EULER, 1:RK_THIRD_ORDER, 2: RK_CLASSIC_FOURTH_ORDER, 3:	BACKWARD_EULER, 4: IMPLICIT_MIDPOINT, 5: CRANK_NICOLSON,\n"
			"6: SDIRK_TWO_STAGES, 7: HEUN_EULER, 8:	BOGACKI_SHAMPINE, 9: DOPRI, 10: FEHLBERG, 11: CASH_KARP,\n"
			"12: RUNGE_KUTTA_5STAGE (NATriuM), 13: THETA_METHOD (NATriuM), 14: EXPONENTIAL"
			<< endl;

	double cFL = 10;
	size_t integrator_id = 0;
	bool automatic_decrease = true;
	if (argc != 6) {
		assert(argc >= 7);
		cFL = atof(argv[6]);
		automatic_decrease = false;
		// automatic crank-nicolson
		integrator_id = 5;
		if (argc > 7){
			assert (argc == 8);
			integrator_id = atoi(argv[7]);
		}
	}

	const size_t cfg = atoi(argv[1]);
	const double Ma = atof(argv[2]);
	const double gamma = atof(argv[3]);
	const size_t refinementLevel = atoi(argv[4]);
	const size_t orderOfFiniteElement = atoi(argv[5]);
	pout << "Ma = " << Ma << ", gamma = " << gamma << endl;
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
	pout << fName.str().c_str() << endl;
	std::ofstream resultFile(fName.str().c_str());
	resultFile
			<< "#<configuration id> <Ma-Number> <Gamma> <Refinement level>  eps     alpha      Lx       h       a        b     u_mean    sigma  tau   Psi_s"
			<< endl;

	pout << "Starting configuration " << cfg << "..." << endl;

	/// get geometry parameters
	double Lx = configurations[cfg][0];
	double h = configurations[cfg][1];
	double a = configurations[cfg][2];
	double b = configurations[cfg][3];
	assert(b < h);
	double sigma = sqrt(0.5 * b * b);
	double epsilon = (h - b) / sigma;
	double alpha = h / Lx;
	double CFL = 0.4;

	/// create CFD problem
	double cell_aspect_ratio = 0.25;
	const double viscosity = u_a * h / Re;
	boost::shared_ptr<ProblemDescription<2> > sinusFlow =
			boost::make_shared<SinusoidalShear2D>(viscosity, u_a, refinementLevel, Lx,
					h, b, cell_aspect_ratio);
	const double dt = CFDSolverUtilities::calculateTimestep<2>(
			*sinusFlow->getMesh(), orderOfFiniteElement, D2Q9(scaling),
			cFL);

	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("OUTPUT_DIR");
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(10);
	configuration->setOutputCheckpointInterval(100000000);
	configuration->setOutputSolutionInterval(1000000);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setCFL(CFL);
	double tau = viscosity / (dt * (scaling * scaling) / 3.0);

	// automatic decrease of timestep to get tau = 0.5
	/*if ((tau < 0.5) and (automatic_decrease)) {
		double new_dt = viscosity / (scaling * scaling / 3.0);
		tau = viscosity / (new_dt * (scaling * scaling) / 3.0);+
		configuration->setCFL(CFL);
		pout << "Config " << cfg << ": tau too small ("
				<< viscosity / (dt * (scaling * scaling) / 3.0)
				<< "). Automatic decrease of time step size (now CFL = "
				<< cFL * new_dt / dt << " , tau = " << tau << ")." << endl;

	}*/

	if (integrator_id < 12){
		configuration->setTimeIntegrator(OTHER);
		configuration->setDealIntegrator(static_cast<DealIntegratorName>(integrator_id));
	} else {
		assert (integrator_id < 15);
		configuration->setTimeIntegrator(static_cast<TimeIntegratorName>(integrator_id - 12));
	}

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
	while (err > 1e-6) {
		for (size_t i = 0; i < 100; i++) {
			solver.stream();
			solver.collide();
			solver.output(solver.getIteration());
			solver.setIteration(solver.getIteration()+1);
		}
		const distributed_vector & ux = solver.getVelocity().at(0);
		qx = u_a * h
				- PhysicalProperties<2>::meanVelocityX(ux,
						solver.getAdvectionOperator()) * h;
		Psi_s_old = Psi_s;
		Psi_s = 2 * qx / (sigma * u_a) - h / sigma;
		err = fabs(Psi_s - Psi_s_old);
	}
	resultFile << cfg << "  " << Ma << "  " << refinementLevel << "  " << gamma << "  "  << epsilon << "   " << alpha
			<< "   " << Lx << "   " << h << "   " << a << "   " << b << "   "
			<< qx / h << "   " << sigma << "   " << tau << "  " << Psi_s
			<< endl;
	pout << "Flow factor " << Psi_s << endl;
	pout << "Analysis finished." << endl;

	return 0;
}
