/*
 * BGKMultistep_test.cpp
 *
 *  Created on: 21.06.2016
 *      Author: dominik
 */


#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/BGKMultistep.h"

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

BOOST_AUTO_TEST_SUITE(BGKMultistep_test)
BOOST_AUTO_TEST_CASE(BGKMultistep_collideAll_test) {

	pout << "BGKMultistep_collideAll_D2Q9_test..." << endl;

	// create collision model// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	double dt = 0.1;
	double viscosity = tau*dt*(1./3.);

	BGKMultistep multistep(tau, dt, boost::make_shared<D2Q9>());
	BGKStandard bgk(tau, dt, boost::make_shared<D2Q9>());
	multistep.setViscosity(viscosity);

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
	DistributionFunctions fAfterCollisionmulti(f);
	DistributionFunctions fAfterCollisionbgk(f);

	for (int p=0;p<1500;p++)
	{

	multistep.collideAll(fAfterCollisionmulti, rho, u, dof_handler.locally_owned_dofs(),
			false);
	bgk.collideAll(fAfterCollisionbgk, rho, u, dof_handler.locally_owned_dofs(),
			false);

	cout << fAfterCollisionmulti.at(4)(1) << " multi to f :  " << f.at(4)(1) << endl;
	cout << fAfterCollisionbgk.at(4)(1) << " bgk to f :  " << f.at(4)(1) << endl;




	for (size_t i = 0; i < dof_handler.n_dofs(); i++) {
		double rho_bgk = 0, rho_multistep = 0;
		if (rho.in_local_range(i)) {
			for (int g = 0; g < 9; g++) {
				rho_bgk = rho_bgk + fAfterCollisionbgk.at(g)(i);
				rho_multistep = rho_multistep + fAfterCollisionmulti.at(g)(i);
			}

			BOOST_CHECK_SMALL(rho_bgk - rho_multistep, 1e-5);
		}}
	}

	pout << "done" << endl;
}
}}
