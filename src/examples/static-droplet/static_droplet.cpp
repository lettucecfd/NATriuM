/**
 * @file static_droplet.cpp
 * @short Static droplet simulation
 * @date 16.11.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"
#include "deal.II/base/utilities.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/Droplet2D.h"
#include "natriuM/benchmarks/PoiseuilleFlow2D.h"

#include "natrium/utilities/Info.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// ========================================================================
	// READ COMMAND LINE PARAMETERS
	// ========================================================================
	pout
			<< "Usage: ./static_droplet <G> <tau> <refinement_level=3> <p=4> <integrator_id = 1>"
			<< endl;
	assert(argc >= 2);
	double G = std::atof(argv[1]);
	double tau = std::atof(argv[2]);

	pout << "Static droplet with: " << endl;
	pout << "... G:    " << G << endl;
	pout << "... tau:  " << tau << endl;

	size_t refinement_level = 3;
	if (argc >= 4) {
		refinement_level = std::atoi(argv[3]);
	}
	pout << "... N:    " << refinement_level << endl;

	size_t p = 4;
	if (argc >= 5) {
		p = std::atoi(argv[4]);
	}
	pout << "... p:    " << p << endl;

	size_t integrator_id = 1;
	if (argc >= 6) {
		integrator_id = std::atoi(argv[5]);
	}
	pout << "... Int:  " << integrator_id << endl;

	// get integrator
	TimeIntegratorName time_integrator;
	DealIntegratorName deal_integrator;
	string integrator_name;
	CFDSolverUtilities::get_integrator_by_id(integrator_id, time_integrator,
			deal_integrator, integrator_name);
	pout << "... that is the " << integrator_name << endl;
	pout << "----------------" << endl;

	// ========================================================================
	// MAKE FLOW PROBLEM
	// ========================================================================
	const double stencil_scaling = 1.0;
	const double CFL = .4 ;
	const double R0 = 0.125;
	const double W = 0.05;
	const double rho_l = 1.932442e+00 ;//2.470937e-01 ;//1.9;//3; //1.93;
	const double rho_g = 1.564130e-01 ;//4.297662e-02 ;//1; //0.15;
	const double length = 1. ; // 2.
	const double height = 1. ;


	// not the real value, is overwritten later (just needed for initialization at this point)
	double viscosity = 1e-1 ;
//	double viscosity = 1e-6 ;

//	boost::shared_ptr<ProblemDescription<2> > droplet2D = boost::make_shared<Droplet2D>(viscosity, refinement_level, 1.0, rho_l, rho_g, W, R0);
	boost::shared_ptr<ProblemDescription<2> > droplet2D = boost::make_shared<Droplet2D>(viscosity,
			refinement_level, length, height, rho_l, rho_g, W, R0);
	double delta_t = CFDSolverUtilities::calculateTimestep<2>(
			*(droplet2D->getMesh()), p, D2Q9(stencil_scaling), CFL);
	droplet2D->setViscosity(
			(tau - 0.5) * delta_t
					* D2Q9(stencil_scaling).getSpeedOfSoundSquare());

	// ========================================================================
	// CONFIGURE SOLVER
	// ========================================================================
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setSwitchOutputOff(false);
	configuration->setUserInteraction(true);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setOutputTableInterval(1);//10
	configuration->setOutputSolutionInterval(1); //10
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/static-droplet";
	configuration->setOutputDirectory(dirname.str());
	configuration->setConvergenceThreshold(1e-10);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(stencil_scaling);
	configuration->setTimeStepSize(delta_t);
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);
	configuration->setCollisionScheme(BGK_MULTIPHASE);
	configuration->setBGKPseudopotentialG(G);

	// ========================================================================
	// RUN SOLVER
	// ========================================================================
	CFDSolver<2> solver(configuration, droplet2D);
	solver.run();

	// ========================================================================
	// FINAL OUTPUT
	// ========================================================================
	pout << "Flow converged" << endl;
	// To evaluate the spurious velocities we have to use an average of the pre- and post-
	// collision velocities. These are evaluated in a final iteration step.
	solver.stream();
	distributed_vector u0( solver.getVelocity().at(0) );
	distributed_vector u1( solver.getVelocity().at(1) );
	solver.collide();
	u0 += solver.getVelocity().at(0);
	u1 += solver.getVelocity().at(1);
	u0 *= 0.5;
	u1 *= 0.5;
	// calculate norm
	u0.scale(u0);
	u1.scale(u1);
	u0 += u1;
	double max_u_square = u0.linfty_norm();
	pout << "Max velocity: " << std::sqrt(max_u_square)
			<< " (while your convergence threshold was "
			<< configuration->getConvergenceThreshold() << ")." << endl;
	return 0;
}
