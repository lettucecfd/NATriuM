/*
 * BGKIncompressible_test.cpp
 *
 *  Created on: 19.09.2015
 *      Author: dominik
 */

#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/BGKIncompressible.h"

#include <math.h>
#include <exception>

#include "boost/test/included/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"

using std::exception;

using namespace natrium;

BOOST_AUTO_TEST_SUITE(BGKIncompressible_test)

BOOST_AUTO_TEST_CASE(BGKIncompressible_collideAll_test) {
/*
	cout << "BGKIncompressible_collideAll_test..." << endl;

	// create collision model// create collision model
	shared_ptr<Stencil> dqmodel = make_shared<D2Q9>();
	double tau = 0.9;
	double dt = 0.1;
	BGKIncompressible bgk(tau, 0.1, make_shared<D2Q9>());

	// initialize distributions with arbitrary components
	vector<distributed_vector> f;
	distributed_vector rho(10);
	vector<distributed_vector> u;
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		distributed_vector f_i(10);
		for (size_t j = 0; j < 10; j++) {
			f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
					+ pow((0.5 * cos(j)), 2);
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
	bgk.collideAll(fAfterCollision, rho, u);
	for (size_t i = 0; i < 10; i++) {
		vector<double> localF(dqmodel->getQ());
		for (size_t j = 0; j < dqmodel->getQ(); j++) {
			localF.at(j) = f.at(j)(i);
		}
		cout << "Lokales f(1) =  " << localF.at(1) << " Lokales f_1(1) = " << fAfterCollision.at(1)(i) << endl;
		bgk.collideSinglePoint(localF);
		for (size_t j = 0; j < dqmodel->getQ(); j++) {
			//cout << i << " " << j << endl;
			BOOST_CHECK(fabs(localF.at(j) - fAfterCollision.at(j)(i)) < 1e-12);
		}
	}*/

	cout << "done" << endl;
}
BOOST_AUTO_TEST_SUITE_END()

