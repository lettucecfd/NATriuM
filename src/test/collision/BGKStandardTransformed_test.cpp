/**
 * @file BGKStandardTransformed_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/BGKStandardTransformed.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"

#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"

using std::exception;

using namespace natrium;

BOOST_AUTO_TEST_SUITE(BGKStandardTransformed_test)

BOOST_AUTO_TEST_CASE(BGKStandardTransformedConstruction_test) {
	pout << "BGKStandardTransformedConstruction_test..." << endl;

	// create Boltzmann model and set relaxation parameter
	double tau = 0.9;
	double dt = 0.1;

	BOOST_CHECK_NO_THROW(
			BGKStandardTransformed bgkCollision(tau, dt, boost::make_shared<D2Q9>()));

	pout << "done" << endl;
} //BGKStandardTransformedConstruction_test

BOOST_AUTO_TEST_CASE(BGKStandardTransformedSetTimeStep_test) {
	pout << "BGKStandardTransformedSetTimeStep_test..." << endl;

	// create Boltzmann model and set relaxation parameter
	double tau = 0.9;
	double dt = 0.1;
	BGKStandardTransformed bgkCollision(tau, dt, boost::make_shared<D2Q9>());

	// check if viscosity is untouched (viscosity ~ dt*tau)
	double dt_times_tau = tau * dt;
	bgkCollision.setTimeStep(0.2);
	BOOST_CHECK_CLOSE(dt_times_tau, bgkCollision.getRelaxationParameter() * 0.2,
			1e-10);

	pout << "done" << endl;
} //BGKStandardTransformedSetTimeStep_test

///////////////////////////////////
// CHECK EQUILIBRIUM //////////////
///////////////////////////////////

BOOST_AUTO_TEST_CASE(BGKTransformedMoments_test) {
	pout << "BGKTransformedMoments_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	// Notice: This text can be applied to different stencils (independent of D and Q)

	const double TOLERANCE = 1e-12;
	// TODO rounding errors are probably too big -> more stable implementation of the eq distribution

	// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	BGKStandardTransformed bgk(tau, 0.1, boost::make_shared<D2Q9>());

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel->getD());
	macroscopicVelocity(0) = 0.02;
	macroscopicVelocity(1) = -0.4;
	if (dqmodel->getD() == 3) {
		macroscopicVelocity(2) = 1.13;
	}

	// calculate equilibrium distributions
	vector<double> eqDistributions(dqmodel->getQ());
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		eqDistributions.at(i) = bgk.getEquilibriumDistribution(i,
				macroscopicVelocity, macroscopicDensity);
	}
	// test first order moment (=density)
	double moment1 = bgk.calculateDensity(eqDistributions);
	BOOST_CHECK_CLOSE(moment1, macroscopicDensity, TOLERANCE);

	// test second order moments (=impulse)
	numeric_vector moment2(dqmodel->getD());
	BOOST_CHECK(moment2(0) == 0);
	BOOST_CHECK(moment2(1) == 0);
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		moment2 += Math::scalar_vector(eqDistributions.at(i),
				dqmodel->getDirection(i));
	}
	numeric_vector impulse = Math::scalar_vector(macroscopicDensity,
			macroscopicVelocity);
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		BOOST_CHECK_SMALL(moment2(i) - impulse(i), TOLERANCE);
	}

	// test impuls tensor (Haenel, p.183)
	numeric_matrix moment3(dqmodel->getD(), dqmodel->getD());
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		for (size_t j = 0; j < dqmodel->getD(); j++) {
			moment3(i, j) = 0.0;
			for (size_t k = 0; k < dqmodel->getQ(); k++) {
				moment3(i, j) += dqmodel->getDirection(k)(i)
						* dqmodel->getDirection(k)(j) * eqDistributions.at(k);
			}
			if (i == j) {
				moment3(i, j) += dqmodel->getSpeedOfSoundSquare()
						* bgk.getRho0();
			}
		}
	}
	numeric_matrix expectedImpulsTensor(dqmodel->getD(), dqmodel->getD());
	double pressure = moment1 * dqmodel->getSpeedOfSoundSquare();
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		for (size_t j = 0; j < dqmodel->getD(); j++) {
			expectedImpulsTensor(i, j) = macroscopicVelocity(i)
					* macroscopicVelocity(j) * macroscopicDensity;
			if (i == j) {
				expectedImpulsTensor(i, j) += pressure;
			}
		}
	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		for (size_t j = 0; j < dqmodel->getD(); j++) {
			BOOST_CHECK_CLOSE(moment3(i, j), expectedImpulsTensor(i, j),
					TOLERANCE);
		}
	}

	pout << "done" << endl;
} //BGKTransformedMoments_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelAllEqDistributions_test) {
	pout << "D2Q9IncompressibleModelAllEqDistributions_test..." << endl;

	// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	BGKStandardTransformed bgk(tau, 0.1, boost::make_shared<D2Q9>());

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel->getD());
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel->getD() == 3) {
		macroscopicVelocity(2) = 1.13;
	}

	// test equilibrium distributions
	vector<double> feq(dqmodel->getQ());
	bgk.getEquilibriumDistributions(feq, macroscopicVelocity,
			macroscopicDensity);
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		BOOST_CHECK_SMALL(
				feq.at(i)
						- bgk.getEquilibriumDistribution(i, macroscopicVelocity,
								macroscopicDensity), 1e-30);
	}

	pout << "done" << endl;
} //D2Q9IncompressibleModelAllEqDistributions_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelMoments_Scaled_test) {
	pout << "D2Q9IncompressibleModelMoments_Scaled_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	// Notice: This text can be applied to different stencils (independent of D and Q)

	const double TOLERANCE = 1e-12;
	// TODO rounding errors are probably too big -> more stable implementation of the eq distribution

	// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>(5.0);
	double tau = 0.9;
	BGKStandardTransformed bgk(tau, 0.1, dqmodel);

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel->getD());
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel->getD() == 3) {
		macroscopicVelocity(2) = 1.13;
	}

	// calculate equilibrium distributions
	vector<double> eqDistributions(dqmodel->getQ());
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		eqDistributions.at(i) = bgk.getEquilibriumDistribution(i,
				macroscopicVelocity, macroscopicDensity);
	}

	// test first order moment (=density)
	double moment1 = bgk.calculateDensity(eqDistributions);
	BOOST_CHECK_SMALL(moment1 - macroscopicDensity, TOLERANCE);

	// test second order moments (=impulse)
	numeric_vector moment2(dqmodel->getD());
	BOOST_CHECK(moment2(0) == 0);
	BOOST_CHECK(moment2(1) == 0);
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		moment2 += Math::scalar_vector(eqDistributions.at(i),
				dqmodel->getDirection(i));
	}
	numeric_vector impulse = Math::scalar_vector(macroscopicDensity,
			macroscopicVelocity);
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		BOOST_CHECK_SMALL(moment2(i) - impulse(i), TOLERANCE);
	}

	// test impuls tensor (Haenel, p.183)
	numeric_matrix moment3(dqmodel->getD(), dqmodel->getD());
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		for (size_t j = 0; j < dqmodel->getD(); j++) {
			moment3(i, j) = 0.0;
			for (size_t k = 0; k < dqmodel->getQ(); k++) {
				moment3(i, j) += dqmodel->getDirection(k)(i)
						* dqmodel->getDirection(k)(j) * eqDistributions.at(k);
			}
			if (i == j) {
				moment3(i, j) += dqmodel->getSpeedOfSoundSquare()
						* bgk.getRho0();
			}
		}
	}
	numeric_matrix expectedImpulsTensor(dqmodel->getD(), dqmodel->getD());
	double pressure = moment1 * dqmodel->getSpeedOfSoundSquare();
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		for (size_t j = 0; j < dqmodel->getD(); j++) {
			expectedImpulsTensor(i, j) = macroscopicVelocity(i)
					* macroscopicVelocity(j) * macroscopicDensity;
			if (i == j) {
				expectedImpulsTensor(i, j) += pressure;
			}
		}
	}
	for (size_t i = 0; i < dqmodel->getD(); i++) {
		for (size_t j = 0; j < dqmodel->getD(); j++) {
			BOOST_CHECK_CLOSE(moment3(i, j), expectedImpulsTensor(i, j),
					TOLERANCE);
		}
	}

	pout << "done" << endl;
} //D2Q9IncompressibleModelMoments_Scaled_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelAllEqDistributions_Scaled_test) {
	pout << "D2Q9IncompressibleModelAllEqDistributions_Scaled_test..." << endl;

	// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>(5.0);
	double tau = 0.9;
	BGKStandardTransformed bgk(tau, 0.1, dqmodel);

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel->getD());
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel->getD() == 3) {
		macroscopicVelocity(2) = 1.13;
	}

	// test equilibrium distributions
	vector<double> feq(dqmodel->getQ());
	bgk.getEquilibriumDistributions(feq, macroscopicVelocity,
			macroscopicDensity);
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		BOOST_CHECK_SMALL(
				feq.at(i)
						- bgk.getEquilibriumDistribution(i, macroscopicVelocity,
								macroscopicDensity), 1e-30);
	}

	pout << "done" << endl;
} //D2Q9IncompressibleModelAllEqDistributions_Scaled_test

///////////////////////////////////
// CHECK COLLISIONS ///////////////
///////////////////////////////////

BOOST_AUTO_TEST_CASE(BGKStandardTransformedCollisionInvariants_test) {
	pout << "BGKStandardTransformedCollisionInvariants_test..." << endl;

	// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	BGKStandardTransformed bgk(tau, 0.1, dqmodel);

	// initialize distributions with arbitrary components
	vector<double> f(bgk.getQ());
	for (size_t i = 0; i < bgk.getQ(); i++) {
		f.at(i) = 0.5 + sin(1.5 * i) + 0.001 + i / (i + 1);
	}

	// calculate macroscopic entities before and after collision
	double rhoBefore = bgk.calculateDensity(f);
	numeric_vector uBefore = bgk.calculateVelocity(f);
	vector<double> fBefore(f);
	bgk.collideSinglePoint(f);
	double rhoAfter = bgk.calculateDensity(f);
	numeric_vector uAfter = bgk.calculateVelocity(f);

	// Check that f was changed
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		BOOST_CHECK(fabs(fBefore.at(i) - f.at(i)) > 1e-5);
	}
	// Check invariance of density
	BOOST_CHECK_SMALL(rhoBefore - rhoAfter, 1e-15);
	// Check invariance of impulse
	BOOST_CHECK_SMALL(rhoBefore * uBefore(0) - rhoAfter * uAfter(0), 1e-15);
	BOOST_CHECK_SMALL(rhoBefore * uBefore(1) - rhoAfter * uAfter(1), 1e-15);

	// Check collision invariance of equality distribution
	double prescribedDensity = 0.8;
	numeric_vector prescribedVelocity(dqmodel->getD());
	vector<double> feq(dqmodel->getQ());
	bgk.getEquilibriumDistributions(feq, prescribedVelocity, prescribedDensity);
	vector<double> feqAfterCollision(feq);
	bgk.collideSinglePoint(feqAfterCollision);
	for (size_t i = 0; i < dqmodel->getQ(); i++) {
		BOOST_CHECK_SMALL(feq.at(i) - feqAfterCollision.at(i), 1e-15);
	}

	pout << "done" << endl;
} //BGKStandardTransformedCollisionInvariants_test

BOOST_AUTO_TEST_CASE(BGKStandardTransformed_collideAll_test) {

	pout << "BGKStandardTransformed_collideAll_test..." << endl;

	// create collision model// create collision model
	boost::shared_ptr<Stencil> dqmodel = boost::make_shared<D2Q9>();
	double tau = 0.9;
	BGKStandardTransformed bgk(tau, 0.1, boost::make_shared<D2Q9>());

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
	DistributionFunctions fAfterCollision(f);
	bgk.collideAll(fAfterCollision, rho, u, dof_handler.locally_owned_dofs());
	for (size_t i = 0; i < dof_handler.n_dofs(); i++) {
		if (rho.in_local_range(i)) {
			vector<double> localF(dqmodel->getQ());
			for (size_t j = 0; j < dqmodel->getQ(); j++) {
				localF.at(j) = f.at(j)(i);
			}
			bgk.collideSinglePoint(localF);
			for (size_t j = 0; j < dqmodel->getQ(); j++) {
				//pout << i << " " << j << endl;
				BOOST_CHECK_CLOSE(localF.at(j), 0.0 + fAfterCollision.at(j)(i), 1e-10);
			}
		} /* if in local range */
	} /* for all dofs */

	pout << "done." << endl;
} /* BGKStandardTransformed_collideAll_test*/

BOOST_AUTO_TEST_SUITE_END()


