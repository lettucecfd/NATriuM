/**
 * @file D2Q25H.cpp
 * @short 
 * @date 09.10.2018
 * @author dwilde3m, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "D2Q25H.h"

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
const size_t D2Q25H::D = 2;
/// Q
const size_t D2Q25H::Q = 25;

/// constructor
D2Q25H::D2Q25H(double scaling) :
		Stencil(2, 25, makeDirections(scaling), makeWeights(), Stencil_D2Q25H,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
                scaling/sqrt(3)), m_speedOfSoundSquare(
                scaling * scaling/(3)), m_scaling(scaling) {
} //constructor

/// destructor
D2Q25H::~D2Q25H() {
} /// destructor



// make weights
vector<double> D2Q25H::makeWeights() {
    double r = (sqrt(5.)-sqrt(2.))/sqrt(3.);

	double w_0 = (-3-3*r*r*r*r+54*r*r)/(75*r*r);

	double w_m = (9*r*r*r*r-6-27*r*r)/(300*r*r*(r*r-1));
	double w_n = (9-6*r*r*r*r-27*r*r)/(300*(1-r*r));


	double w_0n = w_0*w_n;
	double w_0m = w_0*w_m;
	double w_mm = w_m*w_m;
	double w_mn = w_m*w_n;
	double w_nn = w_n*w_n;
	vector<double> result;
	result += w_0*w_0, w_0m,w_0m, w_0m,w_0m, w_mm, w_mm, w_mm, w_mm, w_0n, w_0n, w_0n, w_0n, w_nn, w_nn, w_nn, w_nn, w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn ;
    return result;
} /// make weights

/// make directions
vector<numeric_vector> D2Q25H::makeDirections(double scaling) {
    const double c_m = scaling*sqrt(5.-sqrt(10.))/sqrt(3.);
    const double c_n = sqrt(5.+sqrt(10.))*scaling/sqrt(3.);
	const double directionsArray[][2] = { { 0.0, 0.0 }, { c_m, 0.0 }, { 0.0,
			c_m }, { -c_m, 0.0 }, { 0.0, -c_m },
			{ c_m, c_m }, { -c_m, c_m }, { -c_m, -c_m },
			{ c_m, -c_m } , { c_n, 0.0 }, { 0.0,
					c_n }, { -c_n, 0.0 }, { 0.0, -c_n },
					{ c_n, c_n }, { -c_n, c_n }, { -c_n, -c_n },
					{ c_n, -c_n }, { c_m, c_n }, { c_m, -c_n }, { -c_m, -c_n }, { -c_m, c_n },
							{ c_n, c_m }, { c_n, -c_m }, { -c_n, -c_m },
							{ -c_n, c_m }};
	vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(2);
		direction(0) = directionsArray[i][0];
		direction(1) = directionsArray[i][1];
		result += direction;
	}
	return result;
} /// make directions

numeric_matrix D2Q25H::makeMomentBasis(vector<numeric_vector> e) {
	numeric_matrix m(Q);
	// from Lallemand and Luo (2000)
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
		m(4, alpha) = (-5 + 3 * ealpha_normsq) * e[alpha](0);
		// jy
		m(5, alpha) = e[alpha](1);
		// qy
		m(6, alpha) = (-5 + 3 * ealpha_normsq) * e[alpha](1);
		// pxx
		m(7, alpha) = e[alpha](0)*e[alpha](0)-e[alpha](1)*e[alpha](1);
		// pxy
		m(8, alpha) = e[alpha](0)*e[alpha](1);
		m(9, alpha) = e[alpha](0)*e[alpha](1);
		m(10, alpha) = e[alpha](0)*e[alpha](1);
		m(11, alpha) = e[alpha](0)*e[alpha](1);
		m(12, alpha) = e[alpha](0)*e[alpha](1);
		m(13, alpha) = e[alpha](0)*e[alpha](1);
		m(14, alpha) = e[alpha](0)*e[alpha](1);
		m(15, alpha) = e[alpha](0)*e[alpha](1);
		m(16, alpha) = e[alpha](0)*e[alpha](1);
		m(17, alpha) = e[alpha](0)*e[alpha](1);
		m(18, alpha) = e[alpha](0)*e[alpha](1);
		m(19, alpha) = e[alpha](0)*e[alpha](1);
		m(20, alpha) = e[alpha](0)*e[alpha](1);
		m(21, alpha) = e[alpha](0)*e[alpha](1);
		m(22, alpha) = e[alpha](0)*e[alpha](1);
		m(23, alpha) = e[alpha](0)*e[alpha](1);
		m(24, alpha) = e[alpha](0)*e[alpha](1);

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
	//cout << "D2Q25H:" << endl;
	//m.print(cout);
	return m;
}

} /* namespace natrium */

