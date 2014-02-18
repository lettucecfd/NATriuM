/**
 * @file D2Q9IncompressibleBGKModel_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "utilities/Math.h"
#include "utilities/BasicNames.h"

using std::exception;

namespace natrium {

BOOST_AUTO_TEST_SUITE(D2Q9IncompressibleModel_test)

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelConstruction_test) {
	cout << "D2Q9IncompressibleConstruction_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	D2Q9IncompressibleModel dqmodel;
	BOOST_CHECK_EQUAL(dqmodel.D, 2);
	BOOST_CHECK_EQUAL(dqmodel.Q, 9);
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSound(), sqrt(1./3.));
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSoundSquare(), 1./3.);

	cout << "done" << endl;
} //D2Q9IncompressibleModelConstruction_test


BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelStatic_test) {
	cout << "D2Q9IncompressibleModelStatic_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	BOOST_CHECK_EQUAL(D2Q9IncompressibleModel::D, 2);
	BOOST_CHECK_EQUAL(D2Q9IncompressibleModel::Q, 9);

	cout << "done" << endl;
} //D2Q9IncompressibleStatic_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelGetter_test) {
	cout << "D2Q9IncompressibleModelGetter_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////

	D2Q9IncompressibleModel dqmodel;

	// Dimensions
	BOOST_CHECK_EQUAL(dqmodel.getD(), 2);
	BOOST_CHECK_EQUAL(dqmodel.getQ(), 9);

	// Weights (getWeight; getWeights)
	BOOST_CHECK_EQUAL(dqmodel.getWeight(0),4./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(1),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(2),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(3),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(4),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(5),1./36.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(6),1./36.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(7),1./36.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(8),1./36.);
	vector<double> weightsArray(9, 1./36.);
	weightsArray.at(0) = 4./9.;
	weightsArray.at(1) = 1./9.;
	weightsArray.at(2) = 1./9.;
	weightsArray.at(3) = 1./9.;
	weightsArray.at(4) = 1./9.;
	const vector<double> gotWeights(dqmodel.getWeights());
	BOOST_CHECK_EQUAL_COLLECTIONS(gotWeights.begin(), gotWeights.end(),
			weightsArray.begin(), weightsArray.end());

	// Directions (getDirection; getDirections)
	// direction0 = (0,0)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(0)(0), 0);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(0)(1), 0);
	// direction1 = (1,0)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(1)(0), 1);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(1)(1), 0);
	// direction2 = (0,1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(2)(0), 0);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(2)(1), 1);
	// direction3 = (-1,0)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(3)(0), -1);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(3)(1), 0);
	// direction4 = (0,-1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(4)(0), 0);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(4)(1), -1);
	// direction5 = (1,1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(5)(0), 1);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(5)(1), 1);
	// direction6 = (-1,1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(6)(0), -1);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(6)(1), 1);
	// direction7 = (-1,-1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(7)(0), -1);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(7)(1), -1);
	// direction8 = (1,-1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(8)(0), 1);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(8)(1), -1);
	const vector<numeric_vector> gotDirection(dqmodel.getDirections());
	// direction0 = (0,0)
	BOOST_CHECK_EQUAL(gotDirection.at(0)(0), 0);
	BOOST_CHECK_EQUAL(gotDirection.at(0)(1), 0);
	// direction1 = (1,0)
	BOOST_CHECK_EQUAL(gotDirection.at(1)(0), 1);
	BOOST_CHECK_EQUAL(gotDirection.at(1)(1), 0);
	// direction2 = (0,1)
	BOOST_CHECK_EQUAL(gotDirection.at(2)(0), 0);
	BOOST_CHECK_EQUAL(gotDirection.at(2)(1), 1);
	// direction3 = (-1,0)
	BOOST_CHECK_EQUAL(gotDirection.at(3)(0), -1);
	BOOST_CHECK_EQUAL(gotDirection.at(3)(1), 0);
	// direction4 = (0,-1)
	BOOST_CHECK_EQUAL(gotDirection.at(4)(0), 0);
	BOOST_CHECK_EQUAL(gotDirection.at(4)(1), -1);
	// direction5 = (1,1)
	BOOST_CHECK_EQUAL(gotDirection.at(5)(0), 1);
	BOOST_CHECK_EQUAL(gotDirection.at(5)(1), 1);
	// direction6 = (-1,1)
	BOOST_CHECK_EQUAL(gotDirection.at(6)(0), -1);
	BOOST_CHECK_EQUAL(gotDirection.at(6)(1), 1);
	// direction7 = (-1,-1)
	BOOST_CHECK_EQUAL(gotDirection.at(7)(0), -1);
	BOOST_CHECK_EQUAL(gotDirection.at(7)(1), -1);
	// direction8 = (1,-1)
	BOOST_CHECK_EQUAL(gotDirection.at(8)(0), 1);
	BOOST_CHECK_EQUAL(gotDirection.at(8)(1), -1);

	cout << "done" << endl;
} //D2Q9IncompressibleGetter_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelMoments_test) {
	cout << "D2Q9IncompressibleModelMoments_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	// Notice: This text can be applied to different stencils (independent of D and Q)


	const double TOLERANCE = 1e-14;
	// TODO rounding errors are probably too big -> more stable implementation of the eq distribution

	D2Q9IncompressibleModel dqmodel;

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel.D);
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel.D == 3){
		macroscopicVelocity(2) = 1.13;
	}

	// calculate equilibrium distributions
	vector<double> eqDistributions(dqmodel.Q);
	for (size_t i = 0; i < dqmodel.Q; i++) {
		eqDistributions.at(i) = dqmodel.getEquilibriumDistribution(i, macroscopicVelocity, macroscopicDensity);
	}

	// test first order moment (=density)
	double moment1 = 0.0;
	for (size_t i = 0; i < dqmodel.Q; i++) {
		moment1 += eqDistributions.at(i);
	}
	BOOST_CHECK_SMALL(moment1 - macroscopicDensity, TOLERANCE);

	// test second order moments (=impulse)
	numeric_vector moment2(dqmodel.D);
	BOOST_CHECK(moment2(0) == 0);
	BOOST_CHECK(moment2(1) == 0);
	for (size_t i = 0; i < dqmodel.Q; i++) {
		moment2 += Math::scalar_vector(eqDistributions.at(i), dqmodel.getDirection(i));
	}
	numeric_vector impulse = Math::scalar_vector(macroscopicDensity, macroscopicVelocity);
	for (size_t i = 0; i < dqmodel.D; i++) {
		BOOST_CHECK_SMALL(moment2(i) - impulse(i), TOLERANCE);
	}

	// test impuls tensor (Haenel, p.183)
	numeric_matrix moment3(dqmodel.D, dqmodel.D);
	for (size_t i = 0; i < dqmodel.D; i++) {
		for (size_t j = 0; j < dqmodel.D; j++) {
			moment3(i,j) = 0.0;
			for (size_t k = 0; k < dqmodel.Q; k++) {
				moment3(i,j) += dqmodel.getDirection(k)(i) * dqmodel.getDirection(k)(j) * eqDistributions.at(k);
			}
		}
	}
	numeric_matrix expectedImpulsTensor(dqmodel.D, dqmodel.D);
	double pressure = moment1 * dqmodel.getSpeedOfSoundSquare();
	for (size_t i = 0; i < dqmodel.D; i++) {
		for (size_t j = 0; j < dqmodel.D; j++) {
			expectedImpulsTensor(i,j) = macroscopicVelocity(i) * macroscopicVelocity(j) * macroscopicDensity;
			if (i==j){
				expectedImpulsTensor(i,j) += pressure;
			}
		}
	}
	for (size_t i = 0; i < dqmodel.D; i++) {
		for (size_t j = 0; j < dqmodel.D; j++){
			BOOST_CHECK_SMALL(moment3(i,j) - expectedImpulsTensor(i,j), TOLERANCE);
		}
	}

	cout << "done" << endl;
} //D2Q9IncompressibleModelMoments_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelAllEqDistributions_test) {
	cout << "D2Q9IncompressibleModelAllEqDistributions_test..." << endl;

	D2Q9IncompressibleModel dqmodel;

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel.D);
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel.D == 3){
		macroscopicVelocity(2) = 1.13;
	}

	// test equilibrium distributions
	vector<double> feq(dqmodel.Q);
	dqmodel.getEquilibriumDistributions(feq, macroscopicVelocity, macroscopicDensity);
	for (size_t i = 0; i < dqmodel.Q; i++){
		BOOST_CHECK_SMALL(feq.at(i) - dqmodel.getEquilibriumDistribution(i, macroscopicVelocity, macroscopicDensity), 1e-30);
	}

	cout << "done" << endl;
} //D2Q9IncompressibleModelAllEqDistributions_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelMacroscopicEntities_test) {
	cout << "D2Q9IncompressibleModelMacroscopicEntities_test..." << endl;

	D2Q9IncompressibleModel dqmodel;

	// initialize distributions with arbitrary components
	vector<double> f(dqmodel.getQ());
	for (size_t i = 0; i < dqmodel.getQ(); i++){
		f.at(i) = 0.5 + abs(sin(i));
	}

	// calculate macroscopic entities
	double rho = dqmodel.calculateDensity(f);
	numeric_vector u1 = dqmodel.calculateVelocity(f);
	numeric_vector u2(dqmodel.getD());
	dqmodel.calculateVelocity(f, rho, u2);

	// compare two ways of calculating the velocity
	for (size_t i = 0; i < dqmodel.getD(); i++){
		BOOST_CHECK_EQUAL(u1(i), u2(i));
	}

	// re-calculate and compare density
	double density = 0.0;
	for (size_t i = 0; i < dqmodel.getQ(); i++){
		density += f.at(i);
	}
	BOOST_CHECK_EQUAL(density, rho);

	cout << "done" << endl;
} //D2Q9IncompressibleModelMacroscopicEntities_test



//////////////////////////////////
// TESTS FOR THE SCALED VERSION //
//////////////////////////////////
const double SCALING = 0.5;

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelConstruction_Scaled_test) {
	cout << "D2Q9IncompressibleConstruction_Scaled_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	D2Q9IncompressibleModel dqmodel(SCALING);
	BOOST_CHECK_EQUAL(dqmodel.D, 2);
	BOOST_CHECK_EQUAL(dqmodel.Q, 9);
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSound(), SCALING*sqrt(1./3.));
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSoundSquare(), SCALING*SCALING*1./3.);

	cout << "done" << endl;
} //D2Q9IncompressibleModelConstruction_Scaled_test


BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelGetter_Scaled_test) {
	cout << "D2Q9IncompressibleModelGetter_Scaled_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////

	D2Q9IncompressibleModel dqmodel(SCALING);

	// Dimensions
	BOOST_CHECK_EQUAL(dqmodel.getD(), 2);
	BOOST_CHECK_EQUAL(dqmodel.getQ(), 9);

	// Weights (getWeight; getWeights)
	BOOST_CHECK_EQUAL(dqmodel.getWeight(0),4./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(1),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(2),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(3),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(4),1./9.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(5),1./36.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(6),1./36.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(7),1./36.);
	BOOST_CHECK_EQUAL(dqmodel.getWeight(8),1./36.);
	vector<double> weightsArray(9, 1./36.);
	weightsArray.at(0) = 4./9.;
	weightsArray.at(1) = 1./9.;
	weightsArray.at(2) = 1./9.;
	weightsArray.at(3) = 1./9.;
	weightsArray.at(4) = 1./9.;
	const vector<double> gotWeights(dqmodel.getWeights());
	BOOST_CHECK_EQUAL_COLLECTIONS(gotWeights.begin(), gotWeights.end(),
			weightsArray.begin(), weightsArray.end());

	// Directions (getDirection; getDirections)
	// direction0 = (0,0)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(0)(0), 0);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(0)(1), 0);
	// direction1 = (1,0)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(1)(0), SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(1)(1), 0);
	// direction2 = (0,1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(2)(0), 0);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(2)(1), SCALING);
	// direction3 = (-1,0)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(3)(0), -SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(3)(1), 0);
	// direction4 = (0,-1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(4)(0), 0);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(4)(1), -SCALING);
	// direction5 = (1,1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(5)(0), SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(5)(1), SCALING);
	// direction6 = (-1,1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(6)(0), -SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(6)(1), SCALING);
	// direction7 = (-1,-1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(7)(0), -SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(7)(1), -SCALING);
	// direction8 = (1,-1)
	BOOST_CHECK_EQUAL(dqmodel.getDirection(8)(0), SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getDirection(8)(1), -SCALING);
	const vector<numeric_vector> gotDirection(dqmodel.getDirections());
	// direction0 = (0,0)
	BOOST_CHECK_EQUAL(gotDirection.at(0)(0), 0);
	BOOST_CHECK_EQUAL(gotDirection.at(0)(1), 0);
	// direction1 = (1,0)
	BOOST_CHECK_EQUAL(gotDirection.at(1)(0), SCALING);
	BOOST_CHECK_EQUAL(gotDirection.at(1)(1), 0);
	// direction2 = (0,1)
	BOOST_CHECK_EQUAL(gotDirection.at(2)(0), 0);
	BOOST_CHECK_EQUAL(gotDirection.at(2)(1), SCALING);
	// direction3 = (-1,0)
	BOOST_CHECK_EQUAL(gotDirection.at(3)(0), -SCALING);
	BOOST_CHECK_EQUAL(gotDirection.at(3)(1), 0);
	// direction4 = (0,-1)
	BOOST_CHECK_EQUAL(gotDirection.at(4)(0), 0);
	BOOST_CHECK_EQUAL(gotDirection.at(4)(1), -SCALING);
	// direction5 = (1,1)
	BOOST_CHECK_EQUAL(gotDirection.at(5)(0), SCALING);
	BOOST_CHECK_EQUAL(gotDirection.at(5)(1), SCALING);
	// direction6 = (-1,1)
	BOOST_CHECK_EQUAL(gotDirection.at(6)(0), -SCALING);
	BOOST_CHECK_EQUAL(gotDirection.at(6)(1), SCALING);
	// direction7 = (-1,-1)
	BOOST_CHECK_EQUAL(gotDirection.at(7)(0), -SCALING);
	BOOST_CHECK_EQUAL(gotDirection.at(7)(1), -SCALING);
	// direction8 = (1,-1)
	BOOST_CHECK_EQUAL(gotDirection.at(8)(0), SCALING);
	BOOST_CHECK_EQUAL(gotDirection.at(8)(1), -SCALING);

	cout << "done" << endl;
} //D2Q9IncompressibleGetter_Scaled_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelMoments_Scaled_test) {
	cout << "D2Q9IncompressibleModelMoments_Scaled_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	// Notice: This text can be applied to different stencils (independent of D and Q)


	const double TOLERANCE = 1e-14;
	// TODO rounding errors are probably too big -> more stable implementation of the eq distribution

	D2Q9IncompressibleModel dqmodel(SCALING);

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel.D);
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel.D == 3){
		macroscopicVelocity(2) = 1.13;
	}

	// calculate equilibrium distributions
	vector<double> eqDistributions(dqmodel.Q);
	for (size_t i = 0; i < dqmodel.Q; i++) {
		eqDistributions.at(i) = dqmodel.getEquilibriumDistribution(i, macroscopicVelocity, macroscopicDensity);
	}

	// test first order moment (=density)
	double moment1 = 0.0;
	for (size_t i = 0; i < dqmodel.Q; i++) {
		moment1 += eqDistributions.at(i);
	}
	BOOST_CHECK_SMALL(moment1 - macroscopicDensity, TOLERANCE);

	// test second order moments (=impulse)
	numeric_vector moment2(dqmodel.D);
	BOOST_CHECK(moment2(0) == 0);
	BOOST_CHECK(moment2(1) == 0);
	for (size_t i = 0; i < dqmodel.Q; i++) {
		moment2 += Math::scalar_vector(eqDistributions.at(i), dqmodel.getDirection(i));
	}
	numeric_vector impulse = Math::scalar_vector(macroscopicDensity, macroscopicVelocity);
	for (size_t i = 0; i < dqmodel.D; i++) {
		BOOST_CHECK_SMALL(moment2(i) - impulse(i), TOLERANCE);
	}

	// test impuls tensor (Haenel, p.183)
	numeric_matrix moment3(dqmodel.D, dqmodel.D);
	for (size_t i = 0; i < dqmodel.D; i++) {
		for (size_t j = 0; j < dqmodel.D; j++) {
			moment3(i,j) = 0.0;
			for (size_t k = 0; k < dqmodel.Q; k++) {
				moment3(i,j) += dqmodel.getDirection(k)(i) * dqmodel.getDirection(k)(j) * eqDistributions.at(k);
			}
		}
	}
	numeric_matrix expectedImpulsTensor(dqmodel.D, dqmodel.D);
	double pressure = moment1 * dqmodel.getSpeedOfSoundSquare();
	for (size_t i = 0; i < dqmodel.D; i++) {
		for (size_t j = 0; j < dqmodel.D; j++) {
			expectedImpulsTensor(i,j) = macroscopicVelocity(i) * macroscopicVelocity(j) * macroscopicDensity;
			if (i==j){
				expectedImpulsTensor(i,j) += pressure;
			}
		}
	}
	for (size_t i = 0; i < dqmodel.D; i++) {
		for (size_t j = 0; j < dqmodel.D; j++){
			BOOST_CHECK_SMALL(moment3(i,j) - expectedImpulsTensor(i,j), TOLERANCE);
		}
	}

	cout << "done" << endl;
} //D2Q9IncompressibleModelMoments_Scaled_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelAllEqDistributions_Scaled_test) {
	cout << "D2Q9IncompressibleModelAllEqDistributions_Scaled_test..." << endl;

	D2Q9IncompressibleModel dqmodel(SCALING);

	// Define macroscopic entities
	double macroscopicDensity = 1.45;
	numeric_vector macroscopicVelocity(dqmodel.D);
	macroscopicVelocity(0) = 2.3;
	macroscopicVelocity(1) = -1.14;
	if (dqmodel.D == 3){
		macroscopicVelocity(2) = 1.13;
	}

	// test equilibrium distributions
	vector<double> feq(dqmodel.Q);
	dqmodel.getEquilibriumDistributions(feq, macroscopicVelocity, macroscopicDensity);
	for (size_t i = 0; i < dqmodel.Q; i++){
		BOOST_CHECK_SMALL(feq.at(i) - dqmodel.getEquilibriumDistribution(i, macroscopicVelocity, macroscopicDensity), 1e-30);
	}

	cout << "done" << endl;
} //D2Q9IncompressibleModelAllEqDistributions_Scaled_test

BOOST_AUTO_TEST_CASE(D2Q9IncompressibleModelMacroscopicEntities_Scaled_test) {
	cout << "D2Q9IncompressibleModelMacroscopicEntities_Scaled_test..." << endl;

	D2Q9IncompressibleModel dqmodel(SCALING);

	// initialize distributions with arbitrary components
	vector<double> f(dqmodel.getQ());
	for (size_t i = 0; i < dqmodel.getQ(); i++){
		f.at(i) = 0.5 + abs(sin(i));
	}

	// calculate macroscopic entities
	double rho = dqmodel.calculateDensity(f);
	numeric_vector u1 = dqmodel.calculateVelocity(f);
	numeric_vector u2(dqmodel.getD());
	dqmodel.calculateVelocity(f, rho, u2);

	// compare two ways of calculating the velocity
	for (size_t i = 0; i < dqmodel.getD(); i++){
		BOOST_CHECK_EQUAL(u1(i), u2(i));
	}

	// re-calculate and compare density
	double density = 0.0;
	for (size_t i = 0; i < dqmodel.getQ(); i++){
		density += f.at(i);
	}
	BOOST_CHECK_EQUAL(density, rho);

	cout << "done" << endl;
} //D2Q9IncompressibleModelMacroscopicEntities_Scaled_test

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */

