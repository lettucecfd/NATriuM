/**
 * @file BGKTransformed_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "collisionmodels/BGKTransformed.h"

#include "boost/test/unit_test.hpp"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "boltzmannmodels/BoltzmannModel.h"

#include "utilities/BasicNames.h"


namespace natrium {

BOOST_AUTO_TEST_SUITE(BGKTransformed_test)


BOOST_AUTO_TEST_CASE(BGKTransformedConstruction_test) {
	cout << "BGKTransformedConstruction_test..." << endl;

	// create Boltzmann model and set relaxation parameter
	boost::shared_ptr<BoltzmannModel> dqmodel = boost::make_shared<D2Q9IncompressibleModel>();
	double tau = 0.9;

	BOOST_CHECK_NO_THROW(BGKTransformed bgkCollision(tau, dqmodel));

	cout << "done" << endl;
} //BGKTransformedConstruction_test


BOOST_AUTO_TEST_CASE(BGKTransformedGetter_test) {
	cout << "BGKTransformedGetter_test..." << endl;

	// create collision model
	boost::shared_ptr<BoltzmannModel> dqmodel = boost::make_shared<D2Q9IncompressibleModel>();
	double tau = 0.9;
	BGKTransformed bgkCollision(tau, dqmodel);

	// get relaxation parameter
	BOOST_CHECK_EQUAL(bgkCollision.getRelaxationParameter(),tau);

	cout << "done" << endl;
} //BGKTransformedGetter_test


BOOST_AUTO_TEST_CASE(BGKTransformedCollisionInvariants_test) {
	cout << "BGKTransformedCollisionInvariants_test..." << endl;

	// create collision model
	boost::shared_ptr<BoltzmannModel> dqmodel = boost::make_shared<D2Q9IncompressibleModel>();
	double tau = 0.9;
	BGKTransformed bgkCollision(tau, dqmodel);

	// initialize distributions with arbitrary components
	vector<double> f(dqmodel->getQ());
	for (size_t i = 0; i < dqmodel->getQ(); i++){
		f.at(i) = 1.5 + sin(1.5*i)+0.001+i/(i+1);
	}

	// calculate macroscopic entities before and after collision
	double rhoBefore = dqmodel->calculateDensity(f);
	numeric_vector uBefore = dqmodel->calculateVelocity(f);
	vector<double> fBefore(f);
	bgkCollision.collide(f);
	double rhoAfter = dqmodel->calculateDensity(f);
	numeric_vector uAfter = dqmodel->calculateVelocity(f);

	// Check that f was changed
	for (size_t i=0; i<dqmodel->getQ(); i++){
		BOOST_CHECK(fabs(fBefore.at(i) - f.at(i)) > 1e-5);
	}
	// Check invariance of density
	BOOST_CHECK_SMALL(rhoBefore-rhoAfter, 1e-15);
	// Check invariance of impulse
	BOOST_CHECK_SMALL(rhoBefore*uBefore(0)- rhoAfter*uAfter(0), 1e-15);
	BOOST_CHECK_SMALL(rhoBefore*uBefore(1)- rhoAfter*uAfter(1), 1e-15);

	// Check collision invariance of equality distribution
	double prescribedDensity = 0.8;
	numeric_vector prescribedVelocity(dqmodel->getD());
	vector<double> feq(dqmodel->getQ());
	dqmodel->getEquilibriumDistributions(feq, prescribedVelocity, prescribedDensity);
	vector<double> feqAfterCollision(feq);
	bgkCollision.collide(feqAfterCollision);
	for (size_t i = 0; i < dqmodel->getQ(); i++){
		BOOST_CHECK_SMALL(feq.at(i) - feqAfterCollision.at(i), 1e-15);
	}

	cout << "done" << endl;
} //BGKTransformedCollisionInvariants_test


BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
