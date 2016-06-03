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

			<< "Usage: ./shear-layer <refinement_level=3> <p=4> <collision-id=0 (BGK: 0, KBC: 1)> <semi-lagrange=0> <integrator-id=1> <CFL=0.4> <stencil_scaling=1.0>"
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

	size_t collision_id = 0;
	if (argc >= 4) {
		collision_id = std::atoi(argv[3]);
	}
	pout << "... Coll:  " << collision_id << endl;

	size_t semi_lagrange = 0;
	if (argc >= 5) {
		semi_lagrange = std::atoi(argv[4]);
	}
	pout << "... Filter:  " << semi_lagrange << endl;

	size_t integrator_id = 1;
	if (argc >= 6) {
		integrator_id = std::atoi(argv[5]);
	}
	pout << "... Int:  " << integrator_id << endl;

	double CFL = .4;
	if (argc >= 7) {
		CFL = std::atof(argv[6]);
	}
	pout << "... CFL:  " << CFL << endl;

	double stencil_scaling = 1.0;
	if (argc >= 8) {
		stencil_scaling = std::atof(argv[7]);
	}
	pout << "... stencil_scaling:  " << stencil_scaling << endl;

	// get integrator
	TimeIntegratorName time_integrator;
	DealIntegratorName deal_integrator;
	string integrator_name;
	CFDSolverUtilities::get_integrator_by_id(integrator_id, time_integrator,
			deal_integrator, integrator_name);

	pout << "... that is the " << integrator_name << endl;
	pout << "-------------------------------------" << endl;

	// ========================================================================
	// MAKE FLOW PROBLEM
	// ========================================================================

	const double u0 = 0.04;
	const double kappa = 80;
	const double Re = 30000;
	double viscosity = u0 * 1.0 / Re;
	const double t_c = 2.0 / u0; //twice the eddy turnover time

	boost::shared_ptr<ProblemDescription<2> > shear_layer = boost::make_shared<
			ShearLayer2D>(viscosity, refinement_level, u0, kappa);
	/*double delta_t = CFDSolverUtilities::calculateTimestep<2>(
			*(shear_layer->getMesh()), p, D2Q9(stencil_scaling), CFL);*/

	// **** Grid properties ****
	pout << "**** Grid properties ****" << endl;
	int noCellsInOneDir	= p * pow( 2, refinement_level + 1 );
	pout << "Mesh resolution: " << noCellsInOneDir << "x" << noCellsInOneDir << endl;
	pout << "Number of grid points: " << pow(noCellsInOneDir, 2) << endl;
	pout << "-------------------------------------" << endl;

	// ========================================================================
	// CONFIGURE SOLVER
	// ========================================================================
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setRestartAtIteration(500);
	configuration->setSwitchOutputOff(false);
	configuration->setUserInteraction(false);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setOutputTableInterval(100);	//10
	configuration->setOutputSolutionInterval(100); //10
	configuration->setOutputCheckpointInterval(100);
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/shear-layer/N" << refinement_level
			<< "-p" << p << "-sl" << semi_lagrange << "-coll" << collision_id << "-int" << integrator_id
			<< "-CFL" << CFL << "-scaling" << stencil_scaling;
	configuration->setOutputDirectory(dirname.str());
	configuration->setConvergenceThreshold(1e-10);
	//configuration->setNumberOfTimeSteps(500000);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(stencil_scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(t_c);
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);
	configuration->setOutputTurbulenceStatistics(true);

	if (collision_id == 1) {
		configuration->setCollisionScheme(KBC_STANDARD);
	}

	if (semi_lagrange == 1)
		configuration->setAdvectionScheme(SEMI_LAGRANGIAN);

	/*} else if (semi_lagrange == 2) {
		configuration->setFiltering(true);
		configuration->setFilteringScheme(NEW_FILTER);
	}*/
	//configuration->setFiltering(true);

	pout << "Simulation end time will be t_c = " << t_c << endl;
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
