/**
 * @file BGKPseudopotential_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/collision/BGK.h"

#include "boost/test/unit_test.hpp"

#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/BGKPseudopotential.h"
#include "natrium/advection/SEDGMinLee.h"

#include "natrium/utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(BGKPseudopotential_test)

BOOST_AUTO_TEST_CASE(BGKPseudopotential_collideAll_test) {

	/*
	cout << "BGKPseudopotential_collideAll_test..." << endl;


	// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>(5.0);
	double tau = 0.9;
	double dt = 0.1;
	BGKPseudopotential bgk(tau, 0.1, dqmodel);

	// create problem
	const size_t orderOfFiniteElement = 2;
	const size_t refinementLevel = 2;
	PeriodicTestDomain2D periodic(refinementLevel);
	// advection operator is assigned later as it has to be created, first
	boost::shared_ptr<SEDGMinLee<2> > sedgMinLee = boost::make_shared<SEDGMinLee<2> >(
			periodic.getMesh(), periodic.getBoundaries(),
			orderOfFiniteElement, dqmodel);

	bgk.setAdvectionOperator(sedgMinLee);


	// initialize distributions with arbitrary components
	vector<distributed_vector> f;
	size_t nof_dofs = sedgMinLee->getNumberOfDoFs();
	distributed_vector rho(nof_dofs);
	vector<distributed_vector> u;
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		distributed_vector f_i(nof_dofs);
		for (size_t j = 0; j < nof_dofs; j++) {
			f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
					+ pow((0.5 * cos(j)), 2);
		}
		f.push_back(f_i);
	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		distributed_vector u_i(nof_dofs);
		for (size_t j = 0; j < nof_dofs; j++) {
			u_i(j) = 0;
		}
		u.push_back(u_i);
	}

	// TODO (KNUT) make good test
	// collide and compare to previous collision function
	DistributionFunctions fAfterCollision(f);
	bgk.collideAll(fAfterCollision, rho, u);
	for (size_t i = 0; i < nof_dofs; i++) {
		vector<double> localF(dqmodel->getQ());
		for (size_t j = 0; j < dqmodel->getQ(); j++) {
			localF.at(j) = f.at(j)(i);
		}
		bgk.collideSinglePoint(localF);
		for (size_t j = 0; j < dqmodel->getQ(); j++) {
			//cout << i << " " << j << endl;
			BOOST_CHECK(fabs(localF.at(j) - fAfterCollision.at(j)(i)) < 1e-10);
		}
	}

	cout << "done." << endl; */
} /* BGKPseudopotential_collideAll_test*/

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
