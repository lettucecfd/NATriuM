/*
 * step-6.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: bajat
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/Stencil.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/benchmarks/TaylorGreenVortex3D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-16 ..." << endl;

	/////////////////////////////////////////////////
	// read from command line
	//////////////////////////////////////////////////

	pout
			<< "Usage: ./step-16 <Re=800> <refinement_level=3> <p=4> <collision-id=0 (BGK: 0, KBC: 1)> "
					"<filter=0 (no: 0, exp: 1, new: 2> <integrator-id=1> <CFL=1.0> <stencil-id=0 (D3Q15: 0; D3Q19: 1; D3Q27: 2)>"
			<< endl;

	double Re = 800;
	if (argc >= 2) {
		Re = std::atof(argv[1]);
	}
	pout << "... Re:  " << Re << endl;

	size_t refinement_level = 3;
	if (argc >= 3) {
		refinement_level = std::atoi(argv[2]);
	}
	pout << "... N:    " << refinement_level << endl;

	size_t p = 4;
	if (argc >= 4) {
		p = std::atoi(argv[3]);
	}
	pout << "... p:    " << p << endl;

	size_t collision_id = 0;
	if (argc >= 5) {
		collision_id = std::atoi(argv[4]);
	}
	pout << "... Coll:  " << collision_id << endl;

	size_t filter_id = 0;
	if (argc >= 6) {
		filter_id = std::atoi(argv[5]);
	}
	pout << "... Filter:  " << filter_id << endl;

	size_t integrator_id = 1;
	if (argc >= 7) {
		integrator_id = std::atoi(argv[6]);
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

	double CFL = 1.0;
	if (argc >= 8) {
		CFL = std::atof(argv[7]);
	}
	pout << "... CFL:  " << CFL << endl;

	size_t stencil_id = 0;
	if (argc >= 9) {
		stencil_id = std::atoi(argv[8]);
	}
	pout << "... Sten:  " << stencil_id << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	const double U = 2 * M_PI;
	const double viscosity = 1. / Re;
	// C-E-approach: constant stencil scaling
	// specify Mach number
	const double Ma = 0.1;
	const double cs = U / Ma;
	// zunaechst: fixed order of FE

	// chose scaling so that the right Ma-number is achieved
	const double scaling = sqrt(3) * U / Ma;
	const bool init_rho_analytically = true;

	boost::shared_ptr<ProblemDescription<3> > taylorGreen = boost::make_shared<
			TaylorGreenVortex3D>(viscosity, refinement_level, cs, init_rho_analytically);

	// the scaling has to be orders of magnitude greater than the boundary velocity
	boost::shared_ptr<Stencil> st;
	if (stencil_id == 0){
		st = boost::make_shared<D3Q15>(scaling);
	} else if (stencil_id == 1){
		st = boost::make_shared<D3Q19>(scaling);
	} else if (stencil_id == 2){
		st = boost::make_shared<D3Q27>(scaling);
	}
	double dt = CFDSolverUtilities::calculateTimestep<3>(
			*taylorGreen->getMesh(), p, *st, CFL);

	// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/step-TGV3D/Re" << Re << "-ref"
			<< refinement_level << "-p" << p << "-coll" << collision_id << "-f"
			<< filter_id << "-int" << integrator_id << "-CFL" << CFL << "-sten" << stencil_id;
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	//configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(10);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(1000);
	configuration->setSimulationEndTime(10.0);
	configuration->setInitializationScheme(EQUILIBRIUM);
	configuration->setSedgOrderOfFiniteElement(p);
	configuration->setStencilScaling(scaling);
	configuration->setStencil(Stencil_D3Q15);
	if (stencil_id == 1){
		configuration->setStencil(Stencil_D3Q19);
	} else if (stencil_id == 2){
		configuration->setStencil(Stencil_D3Q27);
	}
	//configuration->setCommandLineVerbosity(BASIC);
	configuration->setCFL(CFL);
	if (dt > 0.1) {
		pout << "Timestep too big." << endl;
	}
	if (collision_id == 1) {
		configuration->setCollisionScheme(KBC_STANDARD);
	}
	configuration->setTimeIntegrator(time_integrator);
	configuration->setDealIntegrator(deal_integrator);
	if (filter_id == 1) {
		configuration->setFiltering(EXPONENTIAL_FILTER);
	}

	//configuration->setNumberOfTimeSteps(1.0 / dt);

	CFDSolver<3> solver(configuration, taylorGreen);

	solver.run();

	pout << "step-16 terminated." << endl;

	return 0;

}
