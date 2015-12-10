/*
 * KBCStandard_test.cpp
 *
 *  Created on: 01.10.2015
 *      Author: dominik
 */

#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/KBCStandard.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

using std::exception;

namespace natrium {

BOOST_AUTO_TEST_SUITE(KBCStandard_test)
BOOST_AUTO_TEST_CASE(KBCStandard_collideAll_test) {

	cout << "KBCStandard_collideAll_test..." << endl;

	// create collision model// create collision model
	shared_ptr<Stencil> dqmodel = make_shared<D2Q9>();
	double tau = 0.9;
	double dt = 0.1;
	vector<distributed_vector> f;
	distributed_vector rho(10);
	vector<distributed_vector> u;

	KBCStandard kbc(tau, dt, make_shared<D2Q9>());
	BGKStandard bgk(tau, dt, make_shared<D2Q9>());

	// initialize distributions with arbitrary components
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		distributed_vector f_i(10);
		for (size_t j = 0; j < 10; j++) {
			f_i(j) = i + 1;
		}
		f.push_back(f_i);
	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		distributed_vector u_i(10);
		for (size_t j = 0; j < 10; j++) {
			u_i(j) = 0;
		}
		u.push_back(u_i);
	}

	// collide and compare to previous collision function
	DistributionFunctions fAfterCollisionkbc(f);
	DistributionFunctions fAfterCollisionbgk(f);
	kbc.collideAll(fAfterCollisionkbc, rho, u, false);
	bgk.collideAll(fAfterCollisionbgk, rho, u, false);

	double rho_bgk=0,rho_kbc=0;

	for (int g=0;g<9;g++)
	{
		rho_bgk = rho_bgk + fAfterCollisionbgk.at(g)(0);
		rho_kbc = rho_kbc + fAfterCollisionbgk.at(g)(0);
	}

	BOOST_CHECK_SMALL(rho_bgk-rho_kbc, 1e-2);

	kbc.setTimeStep(0.2);



	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(KBCStandard_TGV_test) {
	const double Re = 800 * atan(1);
	const double Ma = 0.20;

	const double L = 8 * atan(1); // = 2 pi
	const double U = 1;
	const double tmax = 1;

	// scaling of particle velocities
	double scaling = sqrt(3) * U / Ma;
	// Viscosity
	const double viscosity = U * L / Re;
	// starting time
	//const double t0 = 30.0;
	const double t0 = 1.0; // analytic solution won't converge for t0 = 0.0 and adaptive timesteps

	size_t refinementLevel = 3;
	size_t orderOfFiniteElement = 5;

	shared_ptr<TaylorGreenVortex2D> tgv2D =
			make_shared<TaylorGreenVortex2D>(viscosity, refinementLevel,
					1. / sqrt(3.) / Ma);
	shared_ptr<Benchmark<2> > benchmark = tgv2D;

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();

	configuration->setCollisionScheme(KBC_STANDARD);

	double dt = CFDSolverUtilities::calculateTimestep<2>(
			*benchmark->getTriangulation(), orderOfFiniteElement,
			D2Q9(scaling), 0.2);

	configuration->setTimeIntegrator(RUNGE_KUTTA_5STAGE);

	configuration->setTimeStepSize(dt);
	configuration->setSimulationEndTime(tmax);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setSwitchOutputOff(true);

	BenchmarkCFDSolver<2> solverKBC(configuration, benchmark);

	configuration->setCollisionScheme(BGK_STANDARD);

	BenchmarkCFDSolver<2> solverBGK(configuration, benchmark);

	try{

	solverBGK.run();
	solverBGK.getErrorStats()->update();



	solverKBC.run();
	solverKBC.getErrorStats()->update();

	BOOST_CHECK_SMALL(solverBGK.getErrorStats()->getMaxVelocityError() - solverKBC.getErrorStats()->getMaxVelocityError(), 1e-2);

	cout << "done" << endl;

	} catch (std::exception& e) {
		cout << " Error: " << e.what() << endl;}

}
BOOST_AUTO_TEST_SUITE_END()
}


