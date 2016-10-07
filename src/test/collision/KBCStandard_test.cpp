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

#include "boost/test/included/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"
#include "natrium/collision/BGKStandard.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/benchmarks/PeriodicTestDomain3D.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

using std::exception;

using namespace natrium;

BOOST_AUTO_TEST_SUITE(KBCStandard_test)
BOOST_AUTO_TEST_CASE(KBCStandard_collideAll_test) {

	pout << "KBCStandard_collideAll_D2Q9_test..." << endl;

	// create collision model// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	double dt = 0.1;

	KBCStandard kbc(tau, dt, boost::make_shared<D2Q9>());
	BGKStandard bgk(tau, dt, boost::make_shared<D2Q9>());

	// vectors have to be distributed, because otherwise
	// they are recognized as ghost vectors; and ghost
	// do not support writing on individual elements
	PeriodicTestDomain2D test_domain(3);
	dealii::QGaussLobatto<1> quadrature(2);
	dealii::FE_DGQArbitraryNodes<2> fe(quadrature);
	dealii::DoFHandler<2> dof_handler(*(test_domain.getMesh()));
	dof_handler.distribute_dofs(fe);

	// initialize distributions with arbitrary components
	vector<distributed_vector> f;
	distributed_vector rho;
	rho.reinit((dof_handler.locally_owned_dofs()), MPI_COMM_WORLD);
	//rho.compress(dealii::VectorOperation::add);
	vector<distributed_vector> u;
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		distributed_vector f_i(rho);
		for (size_t j = 0; j < dof_handler.n_dofs(); j++) {
			if (rho.in_local_range(j)) {
				f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
						+ pow((0.5 * cos(j)), 2);
			}
		}
		//f_i.compress(dealii::VectorOperation::add);
		f.push_back(f_i);
	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		distributed_vector u_i(rho);
		for (size_t j = 0; j < 10; j++) {
			u_i(j) = 0;
		}
		//u_i.compress(dealii::VectorOperation::add);
		u.push_back(u_i);
	}

	// collide and compare to previous collision function
	DistributionFunctions fAfterCollisionkbc(f);
	DistributionFunctions fAfterCollisionbgk(f);
	kbc.collideAll(fAfterCollisionkbc, rho, u, dof_handler.locally_owned_dofs(),
			false);
	bgk.collideAll(fAfterCollisionbgk, rho, u, dof_handler.locally_owned_dofs(),
			false);

	for (size_t i = 0; i < dof_handler.n_dofs(); i++) {
		double rho_bgk = 0, rho_kbc = 0;
		if (rho.in_local_range(i)) {
			for (int g = 0; g < 9; g++) {
				rho_bgk = rho_bgk + fAfterCollisionbgk.at(g)(i);
				rho_kbc = rho_kbc + fAfterCollisionkbc.at(g)(i);
			}

			BOOST_CHECK_SMALL(rho_bgk - rho_kbc, 1e-5);
		}
	}

	pout << "done" << endl;
}


BOOST_AUTO_TEST_CASE(KBCStandard_collideAllD3Q15_test) {


	pout << "KBCStandard_collideAll_D3Q15_test..." << endl;

	// create collision model// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D3Q15>();
	double tau = 0.9;
	KBCStandard kbc(tau, 0.1, boost::make_shared<D3Q15>());
	BGKStandard bgk(tau, 0.1, boost::make_shared<D3Q15>());

	// vectors have to be distributed, because otherwise
	// they are recognized as ghost vectors; and ghost
	// do not support writing on individual elements
	PeriodicTestDomain3D test_domain(3);
	dealii::QGaussLobatto<1> quadrature(2);
	dealii::FE_DGQArbitraryNodes<3> fe(quadrature);
	dealii::DoFHandler<3> dof_handler(*(test_domain.getMesh()));
	dof_handler.distribute_dofs(fe);

	// initialize distributions with arbitrary components
	vector<distributed_vector> f;
	distributed_vector rho;
	rho.reinit((dof_handler.locally_owned_dofs()), MPI_COMM_WORLD);
	//rho.compress(dealii::VectorOperation::insert);
	vector<distributed_vector> u;
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		distributed_vector f_i(rho);

		for (size_t j = 0; j < dof_handler.n_dofs(); j++) {
			if (rho.in_local_range(j)) {
				if (rho.in_local_range(j)) {
					f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
							+ pow((0.5 * cos(j)), 2);

			}
			}
		}


		//f_i.compress(dealii::VectorOperation::insert);
		f.push_back(f_i);

	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		distributed_vector u_i(rho);
		for (size_t j = 0; j < 10; j++) {
			u_i(j) = 0;
		}
		//u_i.compress(dealii::VectorOperation::insert);
		u.push_back(u_i);

	}

	// collide and compare to previous collision function
	DistributionFunctions fAfterCollisionkbc(f);
	DistributionFunctions fAfterCollisionbgk(f);

	kbc.collideAll(fAfterCollisionkbc, rho, u, dof_handler.locally_owned_dofs(),
			false);

	bgk.collideAll(fAfterCollisionbgk, rho, u, dof_handler.locally_owned_dofs(),
			false);

	// individual calls for detailed error detection

	for (size_t i = 0; i < dof_handler.n_dofs(); i++) {
		double rho_bgk = 0, rho_kbc = 0;
		if (rho.in_local_range(i)) {
			for (int g = 0; g < 15; g++) {
				rho_bgk = rho_bgk + fAfterCollisionbgk.at(g)(i);
				rho_kbc = rho_kbc + fAfterCollisionkbc.at(g)(i);
			}

			BOOST_CHECK_SMALL(rho_bgk - rho_kbc, 1e-5);
		}
	}
	pout << "done" << endl;
} /* BGKStandard_collideAllD3Q15_test*/

BOOST_AUTO_TEST_SUITE_END()

