/*
 * MRTStandard_test.cpp
 *
 *  Created on: 19.09.2015
 *      Author: dominik
 */

#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/MRTStandard.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"

using std::exception;

namespace natrium {

BOOST_AUTO_TEST_SUITE(MRTStandard_test)

BOOST_AUTO_TEST_CASE(MRTStandard_collideAll_test) {

	cout << "MRTStandard_collideAll_test..." << endl;

	// create collision model// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	double dt = 0.1;

	MRTStandard mrt(tau, dt, boost::make_shared<D2Q9>());
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
	rho.compress(dealii::VectorOperation::add);
	vector<distributed_vector> u;
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		distributed_vector f_i(rho);
		for (size_t j = 0; j < dof_handler.n_dofs(); j++) {
			if (rho.in_local_range(j)) {
				f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
						+ pow((0.5 * cos(j)), 2);
			}
		}
		f_i.compress(dealii::VectorOperation::add);
		f.push_back(f_i);
	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		distributed_vector u_i(rho);
		for (size_t j = 0; j < 10; j++) {
			u_i(j) = 0;
		}
		u_i.compress(dealii::VectorOperation::add);
		u.push_back(u_i);
	}

	// collide and compare to previous collision function
	DistributionFunctions fAfterCollisionmrt(f);
	DistributionFunctions fAfterCollisionbgk(f);
	mrt.collideAll(fAfterCollisionmrt, rho, u, dof_handler.locally_owned_dofs(),
			false);
	bgk.collideAll(fAfterCollisionbgk, rho, u, dof_handler.locally_owned_dofs(),
			false);

	for (size_t i = 0; i < dof_handler.n_dofs(); i++) {
		double rho_bgk = 0, rho_mrt = 0;
		if (rho.in_local_range(i)) {
			for (int g = 0; g < 9; g++) {
				rho_bgk = rho_bgk + fAfterCollisionbgk.at(g)(i);
				rho_mrt = rho_mrt + fAfterCollisionmrt.at(g)(i);
			}

			BOOST_CHECK_SMALL(rho_bgk - rho_mrt, 1e-2);
		}
	}


	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

}
