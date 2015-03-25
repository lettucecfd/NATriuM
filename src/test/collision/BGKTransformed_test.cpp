/**
 * @file BGKTransformed_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "collision/BGKTransformed.h"

#include "boost/test/unit_test.hpp"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "boltzmannmodels/D2Q9PseudopotentialModel.h"
#include "boltzmannmodels/BoltzmannModel.h"

#include "advection/SEDGMinLee.h"
#include "../advection/PeriodicTestDomain2D.h"

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(BGKTransformed_test)

BOOST_AUTO_TEST_CASE(BGKTransformedConstruction_test) {
	cout << "BGKTransformedConstruction_test..." << endl;

	// create Boltzmann model and set relaxation parameter
	D2Q9IncompressibleModel dqmodel;
	double tau = 0.9;
	double dt = 0.1;

	BOOST_CHECK_NO_THROW(BGKTransformed bgkCollision(9,tau, dt));

	cout << "done" << endl;
} //BGKTransformedConstruction_test

BOOST_AUTO_TEST_CASE(BGKTransformedGetter_test) {
	cout << "BGKTransformedGetter_test..." << endl;

	// create collision model
	double tau = 0.9;
	double dt = 0.1;
	BGKTransformed bgkCollision(9, tau, dt);

	cout << "done" << endl;
} //BGKTransformedGetter_test

BOOST_AUTO_TEST_CASE(BGKTransformedSetTimeStep_test) {
	cout << "BGKTransformedSetTimeStep_test..." << endl;

	// create collision model
	double tau = 0.9;
	double dt = 0.1;
	BGKTransformed bgkCollision(9, tau, dt);

	// check if viscosity is untouched (viscosity ~ dt*tau)
	double dt_times_tau = tau*dt;
	bgkCollision.setTimeStep(0.2);
	BOOST_CHECK_CLOSE(dt_times_tau, bgkCollision.getRelaxationParameter()*0.2, 1e-10);

	cout << "done" << endl;
} //BGKTransformedGetter_test

BOOST_AUTO_TEST_CASE(BGKTransformedCollisionInvariants_test) {
	cout << "BGKTransformedCollisionInvariants_test..." << endl;

	// create collision model
	D2Q9IncompressibleModel dqmodel;
	double tau = 0.9;
	double dt = 0.1;
	BGKTransformed bgk(9, tau, 0.1);
	Collision<D2Q9IncompressibleModel, BGKTransformed> bgk_collision(dqmodel,
			bgk);

	// initialize distributions with arbitrary components
	vector<double> f(dqmodel.getQ());
	for (size_t i = 0; i < dqmodel.getQ(); i++) {
		f.at(i) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1);
	}

	// calculate macroscopic entities before and after collision
	double rhoBefore = dqmodel.calculateDensity(f);
	numeric_vector uBefore = dqmodel.calculateVelocity(f);
	vector<double> fBefore(f);
	bgk_collision.collideSinglePoint(f);
	double rhoAfter = dqmodel.calculateDensity(f);
	numeric_vector uAfter = dqmodel.calculateVelocity(f);

	// Check that f was changed
	for (size_t i = 0; i < dqmodel.getQ(); i++) {
		BOOST_CHECK(fabs(fBefore.at(i) - f.at(i)) > 1e-5);
	}
	// Check invariance of density
	BOOST_CHECK_SMALL(rhoBefore - rhoAfter, 1e-15);
	// Check invariance of impulse
	BOOST_CHECK_SMALL(rhoBefore * uBefore(0) - rhoAfter * uAfter(0), 1e-15);
	BOOST_CHECK_SMALL(rhoBefore * uBefore(1) - rhoAfter * uAfter(1), 1e-15);

	// Check collision invariance of equality distribution
	double prescribedDensity = 0.8;
	numeric_vector prescribedVelocity(dqmodel.getD());
	vector<double> feq(dqmodel.getQ());
	dqmodel.getEquilibriumDistributions(feq, prescribedVelocity,
			prescribedDensity);
	vector<double> feqAfterCollision(feq);
	bgk_collision.collideSinglePoint(feqAfterCollision);
	for (size_t i = 0; i < dqmodel.getQ(); i++) {
		BOOST_CHECK_SMALL(feq.at(i) - feqAfterCollision.at(i), 1e-15);
	}

	cout << "done" << endl;
} //BGKTransformedCollisionInvariants_test

BOOST_AUTO_TEST_CASE(BGKTransformed_collideAll_test) {

	cout << "BGKTransformed_collideAll_test..." << endl;

	// create collision model
	D2Q9IncompressibleModel dqmodel;
	double tau = 0.9;
	double dt = 0.1;
	BGKTransformed bgk(9, tau, dt);
	Collision<D2Q9IncompressibleModel, BGKTransformed> bgk_collision(dqmodel,
			bgk);

	// initialize distributions with arbitrary components
	vector<distributed_vector> f;
	distributed_vector rho(10);
	vector<distributed_vector> u;
	for (size_t i = 0; i < dqmodel.getQ(); i++) {
		distributed_vector f_i(10);
		for (size_t j = 0; j < 10; j++) {
			f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
					+ pow((0.5 * cos(j)), 2);
		}
		f.push_back(f_i);
	}
	for (size_t i = 0; i < dqmodel.getD(); i++) {
		distributed_vector u_i(10);
		for (size_t j = 0; j < 10; j++) {
			u_i(j) = 0;
		}
		u.push_back(u_i);
	}

	// collide and compare to previous collision function
	DistributionFunctions fAfterCollision(f);
	bgk_collision.collideAll(fAfterCollision, rho, u);
	for (size_t i = 0; i < 10; i++) {
		vector<double> localF(dqmodel.getQ());
		for (size_t j = 0; j < dqmodel.getQ(); j++) {
			localF.at(j) = f.at(j)(i);
		}
		bgk_collision.collideSinglePoint(localF);
		for (size_t j = 0; j < dqmodel.getQ(); j++) {
			//cout << i << " " << j << endl;
			BOOST_CHECK(fabs(localF.at(j) - fAfterCollision.at(j)(i)) < 1e-15);
		}
	}

	cout << "done." << endl;
} /* BGKTransformed_collideAll_test*/

BOOST_AUTO_TEST_CASE(BGKTransformed_collideAll_PPModel_test) {

	cout << "BGKTransformed_collideAll_PPModel_test..." << endl;

	// create collision model
	const size_t orderOfFiniteElement = 2;
	const double dt = 0.01;
	const size_t refinementLevel = 2;
	PeriodicTestDomain2D periodic(refinementLevel);
	// advection operator is assigned later as it has to be created, first
	shared_ptr<D2Q9PseudopotentialModel> dqmodel = make_shared<
			D2Q9PseudopotentialModel>(1.0, dt);
	shared_ptr<SEDGMinLee<2> > sedgMinLee = make_shared<SEDGMinLee<2> >(
			periodic.getTriangulation(), periodic.getBoundaries(),
			orderOfFiniteElement, dqmodel);
	dqmodel->setAdvectionOperator(sedgMinLee);
	double tau = 0.9;
	BGKTransformed bgk(9, tau, dt);
	Collision<D2Q9PseudopotentialModel, BGKTransformed> bgk_collision(*dqmodel,
			bgk);

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
	bgk_collision.collideAll(fAfterCollision, rho, u);
	for (size_t i = 0; i < nof_dofs; i++) {
		vector<double> localF(dqmodel->getQ());
		for (size_t j = 0; j < dqmodel->getQ(); j++) {
			localF.at(j) = f.at(j)(i);
		}
		bgk_collision.collideSinglePoint(localF);
		for (size_t j = 0; j < dqmodel->getQ(); j++) {
			//cout << i << " " << j << endl;
			BOOST_CHECK(fabs(localF.at(j) - fAfterCollision.at(j)(i)) < 1e-10);
		}
	}

	cout << "done." << endl;
} /* BGKTransformed_collideAll_PPModel_test*/

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
