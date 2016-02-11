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

#include "natrium/benchmarks/ShearLayer2D.h"

#include "natrium/utilities/Info.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// ========================================================================
	// READ COMMAND LINE PARAMETERS
	// ========================================================================
	pout
			<< "Usage: ./shear-layer <refinement_level=3> <p=4> <integrator_id = 1>"
			<< endl;


	size_t refinement_level = 3;
	if (argc >= 2) {
		refinement_level = std::atoi(argv[1]);
	}
	pout << "... N:    " << refinement_level << endl;

	size_t p = 4;
	if (argc >= 3) {
		p = std::atoi(argv[2]);
	}
	pout << "... p:    " << p << endl;

	size_t integrator_id = 1;
	if (argc >= 4) {
		integrator_id = std::atoi(argv[3]);
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
	const double u0 = 0.04;
	const double kappa = 80;
	const double Re = 30000;
	double viscosity = u0 * 1.0 / Re;

	boost::shared_ptr<ProblemDescription<2> > shear_layer = boost::make_shared<ShearLayer2D>(viscosity,
			refinement_level, u0, kappa);
	double delta_t = CFDSolverUtilities::calculateTimestep<2>(
			*(shear_layer->getMesh()), p, D2Q9(stencil_scaling), CFL);
	// ========================================================================
	// CONFIGURE SOLVER
	// ========================================================================
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setSwitchOutputOff(false);
	configuration->setUserInteraction(true);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setOutputTableInterval(10);//10
	configuration->setOutputSolutionInterval(100); //10
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/shear-layer";
	configuration->setOutputDirectory(dirname.str());
	configuration->setConvergenceThreshold(1e-10);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(stencil_scaling);
	configuration->setTimeStepSize(delta_t);
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);

	// ========================================================================
	// RUN SOLVER
	// ========================================================================
	CFDSolver<2> solver(configuration, shear_layer);

	solver.run();

	// ========================================================================
	// FINAL OUTPUT
	// ========================================================================
	pout << "Flow converged" << endl;
	return 0;
}
