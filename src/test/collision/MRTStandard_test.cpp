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
	shared_ptr<Stencil> dqmodel = make_shared<D2Q9>();
	double tau = 0.9;
	double dt = 0.1;
	MRTStandard mrt(tau, 0.1, make_shared<D2Q9>());
	BGKStandard bgk(tau, 0.1, make_shared<D2Q9>());

	cout << " Q= " << dqmodel->getQ() << endl;

	// initialize distributions with arbitrary components
	vector<distributed_vector> f;
	distributed_vector rho(10);
	vector<distributed_vector> u;
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
	mrt.collideAll(fAfterCollisionkbc, rho, u, false);
	bgk.collideAll(fAfterCollisionbgk, rho, u, false);

	double rho_bgk=0,rho_mrt=0;

	for (int g=0;g<9;g++)
	{
		rho_bgk = rho_bgk + fAfterCollisionbgk.at(g)(0);
		rho_mrt = rho_mrt + fAfterCollisionbgk.at(g)(0);
	}

	BOOST_CHECK_SMALL(rho_bgk-rho_mrt, 1e-2);



	cout << "done" << endl;
}
BOOST_AUTO_TEST_SUITE_END()

}
