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

#include "EnstrophySubdomain.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// ========================================================================
	// READ COMMAND LINE PARAMETERS
	// ========================================================================
	pout

			<< "Usage: ./shear-layer <refinement_level=3> <p=4> <collision-id=0 (BGK: 0, KBC: 1)> <semi-lagrange=0> <integrator-id=1> <CFL=0.4> <BDF=0> <filter=0> <filter_s=32> <scaling=1>"
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
	pout << "... Semi-Lagrange:  " << semi_lagrange << endl;

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

	size_t BDF =  0;
	if (argc >= 8) {
		BDF = std::atoi(argv[7]);
	}
	pout << "... BDF: " << BDF << endl;

	size_t filter = 0;
	if (argc >= 9) {
		filter = std::atof(argv[8]);
	}
	pout << "... Filter:  " << filter << endl;

	size_t filter_s = 0;
	if (argc >= 10) {
		filter_s = std::atof(argv[9]);
	}
	pout << "... Filter s:  " << filter_s << endl;

	// scales the speed of sound (scaling=1 <=> Boesch et al.)
	double scaling = 1;
	if (argc >= 11) {
		scaling = std::atof(argv[10]);
	}
	pout << "... stencil scaling:  " << scaling << endl;

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

	// Ma as in KBC paper by Boesch, setup as in Minion and Brown's paper
	const double u0_Boesch = 0.04;
	const double cs_Boesch = 1/sqrt(3);
	const double Ma_Boesch = u0_Boesch/cs_Boesch;
	const double u0 = 1;
	const double cs = u0/Ma_Boesch;
	const double stencil_scaling = sqrt(3) * cs * scaling;
	const double kappa = 80;
	const double viscosity = 0.0001;
	const double t_c = 1;
	const double perturbation = 0.05;

	boost::shared_ptr<ProblemDescription<2> > shear_layer = boost::make_shared<
			ShearLayer2D>(viscosity, refinement_level, u0, kappa, perturbation);
	/*double delta_t = CFDSolverUtilities::calculateTimestep<2>(
	 *(shear_layer->getMesh()), p, D2Q9(stencil_scaling), CFL);*/

	// **** Grid properties ****
	pout << "**** Grid properties ****" << endl;
	int noCellsInOneDir = p * pow(2, refinement_level + 1);
	pout << "Mesh resolution: " << noCellsInOneDir << "x" << noCellsInOneDir
			<< endl;
	pout << "Number of grid points: " << pow(noCellsInOneDir, 2) << endl;
	pout << "-------------------------------------" << endl;

	// ========================================================================
	// CONFIGURE SOLVER
	// ========================================================================
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setRestartAtIteration(500);
	configuration->setSwitchOutputOff(false);
	configuration->setUserInteraction(false);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setOutputTableInterval(10);	//10
	configuration->setOutputSolutionInterval(1000); //10
	configuration->setOutputCheckpointInterval(1e9);
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/shear-layer-MinionBrown/N" << refinement_level
			<< "-p" << p << "-sl" << semi_lagrange << "-coll" << collision_id
			<< "-int" << integrator_id << "-CFL" << CFL << "-BDF"
			<< BDF << "-filter" << filter << "-filt_s" << filter_s << "-scaling" << scaling;
	configuration->setOutputDirectory(dirname.str());
	configuration->setConvergenceThreshold(1e-10);
	//configuration->setNumberOfTimeSteps(500000);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(stencil_scaling);
	configuration->setCFL(CFL);
	configuration->setSimulationEndTime(t_c);
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);
	configuration->setOutputGlobalTurbulenceStatistics(true);

	if (collision_id == 1) {
		configuration->setCollisionScheme(KBC_STANDARD);
	}

	if (semi_lagrange == 1)
		configuration->setAdvectionScheme(SEMI_LAGRANGIAN);

	if (filter == 1) {
		configuration->setFiltering(true);
		configuration->setFilteringScheme(EXPONENTIAL_FILTER);
		configuration->setExponentialFilterAlpha(36);
		configuration->setExponentialFilterS(filter_s);
		configuration->setExponentialFilterNc(4);
	}
	if (BDF == 1)
		configuration->setCollisionScheme(BGK_MULTI_BDF2);

	pout << "Simulation end time will be t_c = " << t_c << endl;
	// ========================================================================
	// RUN SOLVER
	// ========================================================================
	CFDSolver<2> solver(configuration, shear_layer);
	boost::shared_ptr<EnstrophySubdomain> enst = boost::make_shared<EnstrophySubdomain>(solver);
	solver.appendDataProcessor(enst);

	solver.run();

	// ========================================================================
	// FINAL OUTPUT
	// ========================================================================
	pout << "Flow converged" << endl;
	pout << "Enstrophy in subdomain: " << enst->getResult() << endl;
	return 0;
}
