/*
 * Equilibrium_test.cpp
 * */

#include "natrium/collision/BGKStandard.h"

#include <math.h>
#include <array>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/ConfigNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/collision_advanced/Equilibria.h"
#include "natrium/collision_advanced/AuxiliaryCollisionFunctions.h"
#include "natrium/collision_advanced/CollisionSelection.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/solver/CFDSolver.h"

#include "natrium/benchmarks/TaylorGreenVortex3D.h"
#include "../problemdescription/TaylorGreenTest2D.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(CollisionSelection_test_suite)

BOOST_AUTO_TEST_CASE(CollisionSelection2D_test) {

	pout << "CollisionSelection2D_test.." << endl;
	/* prepare everything */
	boost::shared_ptr<ProblemDescription<2> > tgv = boost::make_shared<
			TaylorGreenTest2D>(1.0, 3);
	boost::shared_ptr<SolverConfiguration> cfg = boost::make_shared<
			SolverConfiguration>();
	cfg->setSwitchOutputOff(true);
	CFDSolver<2> solver(cfg, tgv);
	const dealii::IndexSet& dofs =
			solver.getAdvectionOperator()->getLocallyOwnedDofs();
	DistributionFunctions f(solver.getF());
	std::vector<distributed_vector> u;
	distributed_vector rho;
	CFDSolverUtilities::getWriteableVelocity(u, solver.getVelocity(), dofs);
	CFDSolverUtilities::getWriteableDensity(rho, solver.getDensity(), dofs);

	// schemes to test
	cfg->setEquilibriumScheme(BGK_EQUILIBRIUM);
	pout << " - BGK Equilibrium" << endl;
	cfg->setCollisionScheme(BGK_STANDARD);
	pout << " ---- BGK Standard" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	cfg->setCollisionScheme(BGK_REGULARIZED);
	pout << " ---- BGK Regularized" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	cfg->setCollisionScheme(MRT_STANDARD);
	pout << " ---- MRT Standard" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	// -----------------
	cfg->setEquilibriumScheme(INCOMPRESSIBLE_EQUILIBRIUM);
	pout << " - Incompressible Equilibrium" << endl;
	cfg->setCollisionScheme(BGK_STANDARD);
	pout << " ---- BGK Standard" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	cfg->setCollisionScheme(BGK_REGULARIZED);
	pout << " ---- BGK Regularized" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	cfg->setCollisionScheme(MRT_STANDARD);
	pout << " ---- MRT Standard" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	// -----------------
	cfg->setEquilibriumScheme(STEADYSTATE_EQUILIBRIUM);
	pout << " - Steady-state Equilibrium" << endl;
	cfg->setCollisionScheme(BGK_STANDARD);
	pout << " ---- BGK Standard" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	cfg->setCollisionScheme(BGK_REGULARIZED);
	pout << " ---- BGK Regularized" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));
	cfg->setCollisionScheme(MRT_STANDARD);
	pout << " ---- MRT Standard" << endl;
	BOOST_CHECK_NO_THROW(
			selectCollision<2>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
					*solver.getStencil(), false));

	pout << "done." << endl;

}
/* CollisionSelection2D_test */

BOOST_AUTO_TEST_CASE(CollisionSelection3D_test) {

	pout << "CollisionSelection3D_test.." << endl;

	// stencils to test
	std::vector<StencilType> st;
	st.push_back(Stencil_D3Q19);
	st.push_back(Stencil_D3Q15);
	st.push_back(Stencil_D3Q27);

	for (size_t k = 0; k < st.size(); k++) {

		/* prepare everything */
		boost::shared_ptr<ProblemDescription<3> > tgv = boost::make_shared<
				TaylorGreenVortex3D>(1.0, 2);
		boost::shared_ptr<SolverConfiguration> cfg = boost::make_shared<
				SolverConfiguration>();
		cfg->setSwitchOutputOff(true);

		cfg->setStencil(st.at(k));
		CFDSolver<3> solver(cfg, tgv);
		const dealii::IndexSet& dofs =
				solver.getAdvectionOperator()->getLocallyOwnedDofs();
		cout << " - D3Q" << solver.getStencil()->getQ() << endl;
		DistributionFunctions f(solver.getF());
		std::vector<distributed_vector> u;
		distributed_vector rho;
		CFDSolverUtilities::getWriteableVelocity(u, solver.getVelocity(), dofs);
		CFDSolverUtilities::getWriteableDensity(rho, solver.getDensity(), dofs);

		// schemes to test
		cfg->setEquilibriumScheme(BGK_EQUILIBRIUM);
		pout << " ---- BGK Equilibrium" << endl;
		cfg->setCollisionScheme(BGK_STANDARD);
		pout << " ------- BGK Standard" << endl;
		BOOST_CHECK_NO_THROW(
				selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
						*solver.getStencil(), false));
		cfg->setCollisionScheme(BGK_REGULARIZED);
		pout << " ------- BGK Regularized" << endl;
		BOOST_CHECK_NO_THROW(
				selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
						*solver.getStencil(), false));
		if (cfg->getStencil() == Stencil_D3Q19) {
			cfg->setCollisionScheme(MRT_STANDARD);
			cfg->setMRTBasis(DHUMIERES_D3Q19);
			cfg->setMRTRelaxationTimes(RELAX_DHUMIERES_PAPER);
			pout << " ------- MRT Standard" << endl;
			BOOST_CHECK_NO_THROW(
					selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
							*solver.getStencil(), false));
		}

		// ----------
		cfg->setEquilibriumScheme(INCOMPRESSIBLE_EQUILIBRIUM);
		pout << " ---- Incompressible Equilibrium" << endl;
		cfg->setCollisionScheme(BGK_STANDARD);
		pout << " ------- BGK Standard" << endl;
		BOOST_CHECK_NO_THROW(
				selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
						*solver.getStencil(), false));
		cfg->setCollisionScheme(BGK_REGULARIZED);
		pout << " ------- BGK Regularized" << endl;
		BOOST_CHECK_NO_THROW(
				selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
						*solver.getStencil(), false));
		if (cfg->getStencil() == Stencil_D3Q19) {
			cfg->setCollisionScheme(MRT_STANDARD);
			cfg->setMRTBasis(DHUMIERES_D3Q19);
			cfg->setMRTRelaxationTimes(RELAX_DHUMIERES_PAPER);
			pout << " ------- MRT Standard" << endl;
			BOOST_CHECK_NO_THROW(
					selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
							*solver.getStencil(), false));
		}

		// ----------
		cfg->setEquilibriumScheme(STEADYSTATE_EQUILIBRIUM);
		pout << " ---- Steady-state Equilibrium" << endl;
		cfg->setCollisionScheme(BGK_STANDARD);
		pout << " ------- BGK Standard" << endl;
		BOOST_CHECK_NO_THROW(
				selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
						*solver.getStencil(), false));
		cfg->setCollisionScheme(BGK_REGULARIZED);
		pout << " ------- BGK Regularized" << endl;
		BOOST_CHECK_NO_THROW(
				selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
						*solver.getStencil(), false));
		if (cfg->getStencil() == Stencil_D3Q19) {
			cfg->setCollisionScheme(MRT_STANDARD);
			cfg->setMRTBasis(DHUMIERES_D3Q19);
			cfg->setMRTRelaxationTimes(RELAX_DHUMIERES_PAPER);
			pout << " ------- MRT Standard" << endl;
			BOOST_CHECK_NO_THROW(
					selectCollision<3>(*cfg, *tgv, f, rho, u, dofs, 1.0, 0.1,
							*solver.getStencil(), false));
		}


	}
	pout << "done." << endl;

} /* CollisionSelection3D_test */

BOOST_AUTO_TEST_SUITE_END()

