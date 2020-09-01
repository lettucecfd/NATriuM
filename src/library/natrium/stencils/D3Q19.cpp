/*
 * D3Q19.cpp

 *
 *  Created on: Sep 16, 2014
 *      Author: bajat
 */

#include "D3Q19.h"
#include <math.h>

// enable + operator for filling vectors
#include "boost/assign/std/vector.hpp"

// enable + operator for filling vectors
using namespace boost::assign;

namespace natrium {

/////////////////////////////
// ASSIGN STATIC VARIABLES //
/////////////////////////////

/// constructor
D3Q19::D3Q19(double scaling) :
		Stencil(3, 19, makeDirections(scaling), makeWeights(), Stencil_D3Q19,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
				scaling * pow(3, -0.5)), m_speedOfSoundSquare(
				scaling * scaling / 3.), m_scaling(scaling) {
} //constructor

/// destructor
D3Q19::~D3Q19() {
} /// destructor

// make weights
vector<double> D3Q19::makeWeights() {
	vector<double> result;
	result += 1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1.
			/ 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1.
			/ 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.;
	return result;
} /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
vector<numeric_vector> D3Q19::makeDirections(double scaling) {
	const double directionsArray[][3] = { { 0.0, 0.0, 0.0 },
			{ scaling, 0.0, 0.0 }, { 0.0, 0.0, scaling },
			{ -scaling, 0.0, 0.0 }, { 0.0, 0.0, -scaling },
			{ 0.0, -scaling, 0.0 }, { 0.0, scaling, 0.0 }, { scaling, 0.0,
					scaling }, { -scaling, 0.0, scaling }, { -scaling, 0.0,
					-scaling }, { scaling, 0.0, -scaling }, { scaling, -scaling,
					0.0 }, { scaling, scaling, 0.0 },
			{ -scaling, scaling, 0.0 }, { -scaling, -scaling, 0.0 }, { 0.0,
					-scaling, scaling }, { 0.0, scaling, scaling }, { 0.0,
					scaling, -scaling }, { 0.0, -scaling, -scaling } };
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

numeric_matrix D3Q19::makeMomentBasis(vector<numeric_vector> e) {
	// from Conrad et al. (2014): A viscosity adpation method for LB simulations; JCP 276
	/*	double M[19][19]  = { {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
	 {-1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1},
	 {1,-2,-2,-2,-2,-2,-2,1,1,1,1,1,1,1,1,1,1,1,1},
	 {0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0},
	 {0,-2,2,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0},
	 {0,0,0,1,-1,0,0,1,-1,-1,1,0,0,0,0,1,-1,1,-1},
	 {0,0,0,-2,2,0,0,1,-1,-1,1,0,0,0,0,1,-1,1,-1},
	 {0,0,0,0,0,1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1},
	 {0,0,0,0,0,-2,2,0,0,0,0,1,-1,-1,1,1,-1,-1,1},
	 {0,2,2,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-2,-2,-2,-2},
	 {0,-2,-2,1,1,1,1,1,1,1,1,1,1,1,1,-2,-2,-2,-2},
	 {0,0,0,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,0,0,0,0},
	 {0,0,0,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,0,0,0,0},
	 {0,0,0,0,0,0,0,1,1,-1,-1,0,0,0,0,0,0,0,0},
	 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1},
	 {0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1,0,0,0,0},
	 {0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1,0,0,0,0},
	 {0,0,0,0,0,0,0,-1,1,1,-1,0,0,0,0,1,-1,1,-1},
	 {0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,-1,1,1,-1}
	 };*/
	numeric_matrix m(Q);
	// as in d'Humieres original MRT paper
	for (size_t alpha = 0; alpha < Q; alpha++) {
		double ealpha_normsq = e[alpha](0) * e[alpha](0)
				+ e[alpha](1) * e[alpha](1) + e[alpha](2) * e[alpha](2);
		// rho
		m(0, alpha) = 1.0;
		// e
		m(1, alpha) = 19*ealpha_normsq - 30;
		// epsilon
		m(2, alpha) = 0.5*(21*ealpha_normsq*ealpha_normsq -53*ealpha_normsq + 24);
		// jx
		m(3, alpha) = e[alpha](0);
		// qx
		m(4, alpha) = (5*ealpha_normsq-9)*e[alpha](0);

		// jy
		m(5, alpha) = e[alpha](1);
		// qy
		m(6, alpha) = (5*ealpha_normsq-9)*e[alpha](1);
		// jz
		m(7, alpha) = e[alpha](2);
		// qz
		m(8, alpha) = (5*ealpha_normsq-9)*e[alpha](2);
		//3pxx
		m(9, alpha) = 3*e[alpha](0)*e[alpha](0) - ealpha_normsq;
		//3pixx
		m(10, alpha) = (3*ealpha_normsq*ealpha_normsq -5)* (3*e[alpha](0)*e[alpha](0) - ealpha_normsq);
		//pww
		m(11, alpha) = e[alpha](1)*e[alpha](1) - e[alpha](2)*e[alpha](2);
		//piww
		m(12, alpha) = (3*ealpha_normsq*ealpha_normsq -5)* (e[alpha](1)*e[alpha](1) - e[alpha](2)*e[alpha](2));
		//pxy
		m(13, alpha) = e[alpha](0) * e[alpha](1);
		//pyz
		m(14, alpha) = e[alpha](1) * e[alpha](2);
		//pxz
		m(15, alpha) = e[alpha](0) * e[alpha](2);
		//mx
		m(16, alpha) = (e[alpha](1)*e[alpha](1) - e[alpha](2)*e[alpha](2))* e[alpha](0);
		//my
		m(17, alpha) = (e[alpha](2)*e[alpha](2) - e[alpha](0)*e[alpha](0))* e[alpha](1);
		//mz
		m(18, alpha) = (e[alpha](0)*e[alpha](0) - e[alpha](1)*e[alpha](1))* e[alpha](2);
/*
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
		m(10, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1); // 2 1 0
		m(11, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](2); // 2 0 1
		m(12, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](1); // 1 2 0
		//m(13, alpha) = e[alpha](0) * e[alpha](1)
		//		* e[alpha](2); // 1 1 1
		m(13, alpha) = e[alpha](0) * e[alpha](2) * e[alpha](2); // 1 0 2
		m(14, alpha) = e[alpha](1) * e[alpha](1) * e[alpha](2); // 0 2 1
		m(15, alpha) = e[alpha](1) * e[alpha](2) * e[alpha](2); // 0 1 2
		// 0 0 3
		// 0 3 0
		// 3 0 0
		// fourth order moments
		m(16, alpha) = e[alpha](1) * e[alpha](1) * e[alpha](2) * e[alpha](2); // 0 2 2
		m(17, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](2) * e[alpha](2); // 2 0 2
		m(18, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](1); // 2 2 0
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
		 m(26, alpha) = e[alpha](0) * e[alpha](1)
		 * e[alpha](1); // 2 2 2*/
		// revise !
	}
	//cout << "D3Q19:" << endl;
	//m.print(cout);
	return m;

}

}/* namespace natrium */

