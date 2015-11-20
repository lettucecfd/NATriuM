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

	KBCStandard kbc(tau, 0.1, make_shared<D2Q9>());



	cout << "done" << endl;

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
		DistributionFunctions fAfterCollision(f);
		kbc.collideAll(fAfterCollision, rho, u, false);


		cout << "done" << endl;
}
BOOST_AUTO_TEST_SUITE_END()

}






