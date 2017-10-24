/*
 * D3Q13.cpp

 *
 *  Created on: Oct 20, 2017
 *      Author: dwilde
 */

#include "D3Q13.h"
#include <math.h>

// enable + operator for filling vectors
#include "boost/assign/std/vector.hpp"

// enable + operator for filling vectors
using namespace boost::assign;

namespace natrium {

/////////////////////////////
// ASSIGN STATIC VARIABLES //
/////////////////////////////

// assign D and Q
// has to be done outside the class, because function calls are not allowed in initialization of statics
/// D
const size_t D3Q13::D = 3;
/// Q
const size_t D3Q13::Q = 13;

/// constructor
D3Q13::D3Q13(double scaling) :
		Stencil(3, 13, makeDirections(scaling), makeWeights(), Stencil_D3Q13,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
						scaling*0.8506508084
), m_speedOfSoundSquare(
						scaling * scaling * 0.7236067977), m_scaling(scaling) {
	if (scaling > 100) {
		LOG(WARNING)
				<< "The D3Q13 stencil is used with scaling > 100. "
						"This may lead to significant round-off errors."
						"(See UnitTest for Moment matrix of D3Q13 model for details.)"
				<< endl;
	}
} //constructor

/// destructor
D3Q13::~D3Q13() {
} /// destructor

// make weights
vector<double> D3Q13::makeWeights() {
	vector<double> result;
	result += 0.4, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ;
	return result;
} /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
vector<numeric_vector> D3Q13::makeDirections(double scaling) {
	const double phi = scaling* ((1.+ sqrt(5))/2.);
//scaling *= sqrt((5.- sqrt(5))/2.);
	const double directionsArray[][3] = { { 0.0, 0.0, 0.0 },
			{ 0.0, scaling, phi }, { 0.0, -scaling, -phi },
			{ 0.0, -scaling, phi }, { 0.0, scaling, -phi },
			{ scaling, phi, 0.0 }, { -scaling, -phi, 0.0 }, { -scaling, phi, 0.0}, { scaling, -phi, 0.0 }, { phi, 0.0, scaling }, { -phi, 0.0, -scaling }, {
					-phi, 0.0, scaling },
			{phi, 0.0, -scaling}};
	vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(D);
		direction(0) = directionsArray[i][0];
		direction(1) = directionsArray[i][1];
		direction(2) = directionsArray[i][2];
		result += direction;
	}
	return result;
} /// make directions

numeric_matrix D3Q13::makeMomentBasis(vector<numeric_vector> e) {
	numeric_matrix m(Q);
/*	for (size_t alpha = 0; alpha < Q; alpha++) {
		double ealpha_normsq = e[alpha][0] * e[alpha][0]
				+ e[alpha][1] * e[alpha][1] + e[alpha][2] * e[alpha][2];
		// according to d'Humieres; badly conditioned already for |e_alpha| ~= 100
		// rho
		m(0, alpha) = 1.0;
		// e
		m(1, alpha) = ealpha_normsq - 2;
		// epsilon
		m(2, alpha) =
				0.5
						* (15 * ealpha_normsq * ealpha_normsq
								- 55 * ealpha_normsq + 32);
		// jx
		m(3, alpha) = e[alpha](0);
		// qx
		m(4, alpha) = 0.5 * (5 * ealpha_normsq - 13) * e[alpha](0);
		// jy
		m(5, alpha) = e[alpha](1);
		// qy
		m(6, alpha) = 0.5 * (5 * ealpha_normsq - 13) * e[alpha](1);
		// jz
		m(7, alpha) = e[alpha](2);
		// qz
		m(8, alpha) = 0.5 * (5 * ealpha_normsq - 13) * e[alpha](2);
		// 3pxx
		m(9, alpha) = 3 * e[alpha](0) * e[alpha](0) - ealpha_normsq;
		// pww
		m(10, alpha) = e[alpha](1) * e[alpha](1) - e[alpha](2) * e[alpha](2);
		// pxy
		m(11, alpha) = e[alpha](0) * e[alpha](1);
		// pyz
		m(12, alpha) = e[alpha](1) * e[alpha](2);
		// pzx */


		/*// density
		 m(0, alpha) = 1.0;
		 // velocity
		 m(1, alpha) = e[alpha](0);
		 m(2, alpha) = e[alpha](1);
		 m(3, alpha) = e[alpha](2);
		 // second order moments
		 m(4, alpha) = e[alpha](0) * e[alpha](0); // 2 0 0
		 m(5, alpha) = e[alpha](1) * e[alpha](1); // 0 2 0
		 m(6, alpha) = e[alpha](2) * e[alpha](2); // 0 0 2
		 m(7, alpha) = e[alpha](0) * e[alpha](1); // 1 1 0
		 m(8, alpha) = e[alpha](0) * e[alpha](2); // 1 0 1
		 m(9, alpha) = e[alpha](1) * e[alpha](2); // 0 1 1
		 // third order moments
		 //m(10, alpha) = e[alpha](0) * e[alpha](0)
		 //		* e[alpha](1); // 2 1 0
		 //m(11, alpha) = e[alpha](0) * e[alpha](0)
		 //		* e[alpha](2); // 2 0 1
		 //m(12, alpha) = e[alpha](0) * e[alpha](1)
		 //		* e[alpha](1); // 1 2 0
		 m(10, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](2); // 1 1 1
		 //m(14, alpha) = e[alpha](0) * e[alpha](2)
		 //		* e[alpha](2); // 1 0 2
		 //m(15, alpha) = e[alpha](1) * e[alpha](1)
		 //		* e[alpha](2); // 0 2 1
		 //m(16, alpha) = e[alpha](1) * e[alpha](2)
		 //		* e[alpha](2); // 0 1 2
		 m(11, alpha) = e[alpha](2) * e[alpha](2) * e[alpha](2);		// 0 0 3
		 m(12, alpha) = e[alpha](1) * e[alpha](1) * e[alpha](1);		// 0 3 0
		 m(13, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](0);		// 3 0 0
		 // fourth order moments
		 //m(17, alpha) = e[alpha](1) * e[alpha](1)
		 //		* e[alpha](2) * e[alpha](2); // 0 2 2
		 //m(18, alpha) = e[alpha](0) * e[alpha](1)
		 //		* e[alpha](2) * e[alpha](2); // 2 0 2
		 //m(19, alpha) = e[alpha](0) * e[alpha](0)
		 //		* e[alpha](1) * e[alpha](1); // 2 2 0
		 // 0 0 4
		 // 0 4 0
		 // 4 0 0
		 // 0 1 3
		 // 0 3 1
		 // 1 0 3
		 // 1 3 0
		 // 3 0 1
		 // 3 1 0
		 // 1 1 2
		 // 1 2 1
		 // 2 1 1
		 // fifth order moments
		 m(20, alpha) = e[alpha](0) * e[alpha](1)
		 * e[alpha](1) * e[alpha](1); // 2 1 1
		 m(21, alpha) = e[alpha](0) * e[alpha](0)
		 * e[alpha](1) * e[alpha](1); // 1 1 2
		 m(22, alpha) = e[alpha](0) * e[alpha](1)
		 * e[alpha](1) * e[alpha](1); // 1 2 1
		 m(23, alpha) = e[alpha](0) * e[alpha](0)
		 * e[alpha](1) * e[alpha](0); // 1 2 2
		 m(24, alpha) = e[alpha](0) * e[alpha](1)
		 * e[alpha](1) * e[alpha](0); // 2 1 2
		 m(25, alpha) = e[alpha](0) * e[alpha](0)
		 * e[alpha](1) * e[alpha](0); // 2 2 1
		 // ...
		 // sixth order moment
		 m(14, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](1)
		 * e[alpha](2) * e[alpha](2); // 2 2 2
		 // revise !
		 */

	//cout << "D3Q13:" << endl;
	//m.print(cout);
	return m;
}

}/* namespace natrium */

