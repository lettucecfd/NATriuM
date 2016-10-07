/**
 * @file D2Q9_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/stencils/D2Q9.h"

#include <math.h>
#include <exception>

#include "boost/test/included/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"

using std::exception;

using namespace natrium;

BOOST_AUTO_TEST_SUITE(D2Q9_test)

BOOST_AUTO_TEST_CASE(D2Q9Construction_test) {
	pout << "D2Q9Construction_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	D2Q9 dqmodel;
	BOOST_CHECK_EQUAL(dqmodel.D, size_t(2));
	BOOST_CHECK_EQUAL(dqmodel.Q, size_t(9));
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSound(), sqrt(1./3.));
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSoundSquare(), 1./3.);

	pout << "done" << endl;
} //D2Q9Construction_test


BOOST_AUTO_TEST_CASE(D2Q9Static_test) {
	pout << "D2Q9Static_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	BOOST_CHECK_EQUAL(D2Q9::D, size_t(2));
	BOOST_CHECK_EQUAL(D2Q9::Q, size_t(9));

	pout << "done" << endl;
} //D2Q9Static_test

BOOST_AUTO_TEST_CASE(D2Q9Getter_test) {
	pout << "D2Q9Getter_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////

	D2Q9 dqmodel;

	// Dimensions
	BOOST_CHECK_EQUAL(dqmodel.getD(), size_t(2));
	BOOST_CHECK_EQUAL(dqmodel.getQ(), size_t(9));

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

	pout << "done" << endl;
} //D2Q9Getter_test


//////////////////////////////////
// TESTS FOR THE SCALED VERSION //
//////////////////////////////////
const double SCALING = 0.5;

BOOST_AUTO_TEST_CASE(D2Q9Construction_Scaled_test) {
	pout << "D2Q9Construction_Scaled_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////
	D2Q9 dqmodel(SCALING);
	BOOST_CHECK_EQUAL(dqmodel.getD(), size_t(2));
	BOOST_CHECK_EQUAL(dqmodel.getQ(), size_t(9));
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSound(), SCALING*sqrt(1./3.));
	BOOST_CHECK_EQUAL(dqmodel.getSpeedOfSoundSquare(), SCALING*SCALING*1./3.);

	pout << "done" << endl;
} //D2Q9Construction_Scaled_test


BOOST_AUTO_TEST_CASE(D2Q9Getter_Scaled_test) {
	pout << "D2Q9Getter_Scaled_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////

	D2Q9 dqmodel(SCALING);

	// Dimensions
	BOOST_CHECK_EQUAL(dqmodel.getD(), size_t(2));
	BOOST_CHECK_EQUAL(dqmodel.getQ(), size_t(9));

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
	BOOST_CHECK_CLOSE(dqmodel.getDirection(0)(0), 0, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(0)(1), 0, 1e-10);
	// direction1 = (1,0)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(1)(0), SCALING, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(1)(1), 0, 1e-10);
	// direction2 = (0,1)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(2)(0), 0, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(2)(1), SCALING, 1e-10);
	// direction3 = (-1,0)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(3)(0), -SCALING, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(3)(1), 0, 1e-10);
	// direction4 = (0,-1)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(4)(0), 0, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(4)(1), -SCALING, 1e-10);
	// direction5 = (1,1)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(5)(0), SCALING, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(5)(1), SCALING, 1e-10);
	// direction6 = (-1,1)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(6)(0), -SCALING, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(6)(1), SCALING, 1e-10);
	// direction7 = (-1,-1)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(7)(0), -SCALING, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(7)(1), -SCALING, 1e-10);
	// direction8 = (1,-1)
	BOOST_CHECK_CLOSE(dqmodel.getDirection(8)(0), SCALING, 1e-10);
	BOOST_CHECK_CLOSE(dqmodel.getDirection(8)(1), -SCALING, 1e-10);
	const vector<numeric_vector> gotDirection(dqmodel.getDirections());
	// direction0 = (0,0)
	BOOST_CHECK_CLOSE(gotDirection.at(0)(0), 0, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(0)(1), 0, 1e-10);
	// direction1 = (1,0)
	BOOST_CHECK_CLOSE(gotDirection.at(1)(0), SCALING, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(1)(1), 0, 1e-10);
	// direction2 = (0,1)
	BOOST_CHECK_CLOSE(gotDirection.at(2)(0), 0, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(2)(1), SCALING, 1e-10);
	// direction3 = (-1,0)
	BOOST_CHECK_CLOSE(gotDirection.at(3)(0), -SCALING, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(3)(1), 0, 1e-10);
	// direction4 = (0,-1)
	BOOST_CHECK_CLOSE(gotDirection.at(4)(0), 0, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(4)(1), -SCALING, 1e-10);
	// direction5 = (1,1)
	BOOST_CHECK_CLOSE(gotDirection.at(5)(0), SCALING, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(5)(1), SCALING, 1e-10);
	// direction6 = (-1,1)
	BOOST_CHECK_CLOSE(gotDirection.at(6)(0), -SCALING, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(6)(1), SCALING, 1e-10);
	// direction7 = (-1,-1)
	BOOST_CHECK_CLOSE(gotDirection.at(7)(0), -SCALING, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(7)(1), -SCALING, 1e-10);
	// direction8 = (1,-1)
	BOOST_CHECK_CLOSE(gotDirection.at(8)(0), SCALING, 1e-10);
	BOOST_CHECK_CLOSE(gotDirection.at(8)(1), -SCALING, 1e-10);

	pout << "done" << endl;
} //D2Q9Getter_Scaled_test

BOOST_AUTO_TEST_CASE(D2Q9MomentTrafo_test) {
	pout << "D2Q9MomentTrafo_test..." << endl;

	/////////////////
	// SANITY TEST //
	/////////////////

	// standard scaling
	D2Q9 dqmodel(1.0);

	numeric_matrix f_to_m(9);
	numeric_matrix m_to_f(9);
	dqmodel.getMomentBasis(f_to_m);
	dqmodel.getInverseMomentBasis(m_to_f);
	numeric_matrix Ident(9);
	m_to_f.mmult(Ident, f_to_m);
	for (size_t i = 0; i < 9; i++){
		for (size_t j = 0; j< 9; j++){
			if (i != j)
				BOOST_CHECK_SMALL(Ident(i,j), 1e-10);
			else
				BOOST_CHECK_SMALL(Ident(i,j)-1.0, 1e-10);
		}
	}

	D2Q9 dqmodel2(1000.0);

	dqmodel.getMomentBasis(f_to_m);
	dqmodel.getInverseMomentBasis(m_to_f);
	m_to_f.mmult(Ident, f_to_m);
	for (size_t i = 0; i < 9; i++){
		for (size_t j = 0; j< 9; j++){
			if (i != j)
				BOOST_CHECK_SMALL(Ident(i,j), 1e-10);
			else
				BOOST_CHECK_SMALL(Ident(i,j)-1.0, 1e-10);
		}
	}


	pout << "done" << endl;
} //D2Q9MomentTrafo_test


BOOST_AUTO_TEST_SUITE_END()


