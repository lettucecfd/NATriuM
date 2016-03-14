/**
 * @file D2Q9.cpp
 * @short 
 * @date 30.08.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "D2Q9.h"

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
const size_t D2Q9::D = 2;
/// Q
const size_t D2Q9::Q = 9;

/// constructor
D2Q9::D2Q9(double scaling) :
		Stencil(2, 9, makeDirections(scaling), makeWeights(), Stencil_D2Q9,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
				scaling * pow(3, -0.5)), m_speedOfSoundSquare(
				scaling * scaling / 3.), m_scaling(scaling) {
} //constructor

/// destructor
D2Q9::~D2Q9() {
} /// destructor

// make weights
vector<double> D2Q9::makeWeights() {
	vector<double> result;
	result += 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1.
			/ 36., 1. / 36.;
	return result;
} /// make weights

/// make directions
vector<numeric_vector> D2Q9::makeDirections(double scaling) {
	const double directionsArray[][2] = { { 0.0, 0.0 }, { scaling, 0.0 }, { 0.0,
			scaling }, { -scaling, 0.0 }, { 0.0, -scaling },
			{ scaling, scaling }, { -scaling, scaling }, { -scaling, -scaling },
			{ scaling, -scaling } };
	vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(2);
		direction(0) = directionsArray[i][0];
		direction(1) = directionsArray[i][1];
		result += direction;
	}
	return result;
} /// make directions

numeric_matrix D2Q9::makeMomentBasis(vector<numeric_vector> e) {
	numeric_matrix m(Q);
	// from Bachelor thesis Christoph Rettinger (Fluid flow simulatios using the LBM with MRT)
	for (size_t alpha = 0; alpha < Q; alpha++) {
		double ealpha_normsq = e[alpha](0)*e[alpha](0) + e[alpha](1)*e[alpha](1) ;
		// rho
		m(0, alpha) = 1.0;
		// e
		m(1, alpha) = -4 + 3*ealpha_normsq;
		// eps
		m(2, alpha) = 4 - 21./2. * ealpha_normsq + 9./2. * ealpha_normsq * ealpha_normsq;
		// jx
		m(3, alpha) = e[alpha](0);
		// qx
		m(4, alpha) = (-5 + 3 * ealpha_normsq) * e[alpha](0)* e[alpha](0);
		// jy
		m(5, alpha) = e[alpha](1);
		// qy
		m(6, alpha) = (-5 + 3 * ealpha_normsq) * e[alpha](1)* e[alpha](1);
		// pxx
		m(7, alpha) = e[alpha](0)*e[alpha](0)-e[alpha](1)*e[alpha](1);
		// pxy
		m(8, alpha) = e[alpha](0)*e[alpha](1);

		/*// density
		m(0, alpha) = 1.0;
		// velocity
		m(1, alpha) = e[alpha](0);
		m(2, alpha) = e[alpha](1);
		// second order moments
		m(3, alpha) = e[alpha](0) * e[alpha](0);
		m(4, alpha) = e[alpha](1) * e[alpha](1);
		m(5, alpha) = e[alpha](0) * e[alpha](1);
		// third order moments
		m(6, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](1); // 1 2
		m(7, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1); // 2 1
		// 3 0
		// 0 3
		// fourth order moment
		m(8, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](1); // 2 2*/
	}
	cout << "D2Q9:" << endl;
	m.print(cout);
	return m;
}

} /* namespace natrium */

