/*
 * D3Q27.cpp

 *
 *  Created on: Nov 20, 2015
 *      Author: kraemer
 */

#include "D3Q27.h"
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
const size_t D3Q27::D = 3;
/// Q
const size_t D3Q27::Q = 27;

/// constructor
D3Q27::D3Q27(double scaling) :
		Stencil(3, 27, makeDirections(scaling), makeWeights(), Stencil_D3Q27,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
				scaling * pow(3, -0.5)), m_speedOfSoundSquare(
				scaling * scaling / 3.), m_scaling(scaling) {
} //constructor

/// destructor
D3Q27::~D3Q27() {
} /// destructor

// make weights
vector<double> D3Q27::makeWeights() {
	vector<double> result;
	result += 8. / 27.,
			2. / 27.,
			2. / 27.,
			2. / 27.,
			2. / 27.,
			2. / 27.,
			2. / 27.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 54.,
			1. / 216.,
			1. / 216.,
			1. / 216.,
			1. / 216.,
			1. / 216.,
			1. / 216.,
			1. / 216.,
			1. / 216.;
	return result;
} /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
vector<numeric_vector> D3Q27::makeDirections(double scaling) {
	const double directionsArray[][3] =
	{
			{ 0.0, 0.0, 0.0 },
			{ scaling, 0.0, 0.0 },
			{ -scaling, 0.0, 0.0 },
			{ 0.0, scaling, 0.0 },
			{ 0.0, -scaling, 0.0 },
			{ 0.0, 0.0, scaling },
			{ 0.0, 0.0, -scaling },
			{ 0.0, scaling,	scaling },
			{ 0.0, -scaling, -scaling },
			{ 0.0, scaling,	-scaling },
			{ 0.0, -scaling, scaling },
			{ scaling, 0.0,	scaling },
			{ -scaling, 0.0, -scaling },
			{ scaling, 0.0,	-scaling },
			{ -scaling, 0.0, scaling },
			{ scaling, scaling,	0.0 },
			{ -scaling, -scaling, 0.0 },
			{ scaling, -scaling, 0.0 },
			{ -scaling, scaling, 0.0 },
			{ scaling, scaling, scaling },
			{ -scaling, -scaling, -scaling },
			{scaling, scaling, -scaling },
			{ -scaling, -scaling, scaling },
			{ scaling, -scaling, scaling },
			{-scaling, scaling, -scaling },
			{ scaling, -scaling, -scaling },
			{ -scaling, scaling, scaling } };
	vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(D);
		direction(0) = directionsArray[i][0];
		direction(1) = directionsArray[i][1];
		direction(2) = directionsArray[i][2];
		result += direction;
	}
	return result;
}
/// make directions

numeric_matrix D3Q27::makeMomentBasis(vector<numeric_vector> e) {

	numeric_matrix m(Q);
	// as in Suga et al. (2015): A D3Q27 MRT LBM for turbulent flows
	for (size_t alpha = 0; alpha < Q; alpha++) {
		double ealpha_normsq = e[alpha][0] * e[alpha][0]
				+ e[alpha][1] * e[alpha][1] + e[alpha][2] * e[alpha][2];
		// rho
		m(0, alpha) = 1;
		// jx
		m(1, alpha) = e[alpha](0);
		// jy
		m(2, alpha) = e[alpha](1);
		// jz
		m(3, alpha) = e[alpha](2);
		// e
		m(4, alpha) = ealpha_normsq;
		// xx
		m(5, alpha) = 3 * e[alpha](0) * e[alpha](0) - ealpha_normsq;
		// ww
		m(6, alpha) = e[alpha](1) * e[alpha](1) - e[alpha](2) * e[alpha](2);
		// xy
		m(7, alpha) = e[alpha](0) * e[alpha](1);
		// yz
		m(8, alpha) = e[alpha](1) * e[alpha](2);
		// zx
		m(9, alpha) = e[alpha](0) * e[alpha](2);
		// phi_x
		m(10, alpha) = 3 * ealpha_normsq * e[alpha](0);
		// phi_y
		m(11, alpha) = 3 * ealpha_normsq * e[alpha](1);
		// phi_z
		m(12, alpha) = 3 * ealpha_normsq * e[alpha](2);
		// psi_x
		m(13, alpha) = 9. / 2. * ealpha_normsq * ealpha_normsq * e[alpha](0);
		// psi_y
		m(14, alpha) = 9. / 2. * ealpha_normsq * ealpha_normsq * e[alpha](1);
		// psi_z
		m(15, alpha) = 9. / 2. * ealpha_normsq * ealpha_normsq * e[alpha](2);
		// eps
		m(16, alpha) = 3. / 2. * ealpha_normsq * ealpha_normsq;
		// e3
		m(17, alpha) = 9. / 2. * ealpha_normsq * ealpha_normsq * ealpha_normsq;
		// xxe
		m(18, alpha) = m(5, alpha) * ealpha_normsq;
		// wwe
		m(19, alpha) = m(6, alpha) * ealpha_normsq;
		// xye
		m(20, alpha) = m(7, alpha) * ealpha_normsq;
		// yze
		m(21, alpha) = m(8, alpha) * ealpha_normsq;
		// zxe
		m(22, alpha) = m(9, alpha) * ealpha_normsq;
		// tx
		m(23, alpha) = e[alpha](0)
				* (e[alpha](1) * e[alpha](1) - e[alpha](2) * e[alpha](2));
		// ty
		m(24, alpha) = e[alpha](1)
				* (e[alpha](2) * e[alpha](2) - e[alpha](0) * e[alpha](0));
		// tz
		m(25, alpha) = e[alpha](2)
				* (e[alpha](0) * e[alpha](0) - e[alpha](1) * e[alpha](1));
		// xyz
		m(26, alpha) = e[alpha](0)*e[alpha](1)*e[alpha](2);

		/*		// density
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
		 m(10, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1);	// 2 1 0
		 m(11, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](2);	// 2 0 1
		 m(12, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](1);	// 1 2 0
		 m(13, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](2);	// 1 1 1
		 m(14, alpha) = e[alpha](0) * e[alpha](2) * e[alpha](2);	// 1 0 2
		 m(15, alpha) = e[alpha](1) * e[alpha](1) * e[alpha](2);	// 0 2 1
		 m(16, alpha) = e[alpha](1) * e[alpha](2) * e[alpha](2);	// 0 1 2
		 // 0 0 3
		 // 0 3 0
		 // 3 0 0
		 // fourth order moments
		 m(17, alpha) = e[alpha](1) * e[alpha](1) * e[alpha](2) * e[alpha](2);// 0 2 2
		 m(18, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](2) * e[alpha](2);// 2 0 2
		 m(19, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](1);// 2 2 0
		 // 0 0 4
		 // 0 4 0
		 // 4 0 0
		 // 0 1 3
		 // 0 3 1
		 // 1 0 3
		 // 1 3 0
		 // 3 0 1
		 // 3 1 0
		 m(20, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](2);// 2 1 1
		 m(21, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](2) * e[alpha](2);// 1 1 2
		 m(22, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](1) * e[alpha](2);// 1 2 1
		 // fifth order moments
		 m(23, alpha) = e[alpha](0) * e[alpha](1) * e[alpha](1) * e[alpha](2)
		 * e[alpha](2);	// 1 2 2
		 m(24, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](2)
		 * e[alpha](2);	// 2 1 2
		 m(25, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](1)
		 * e[alpha](2);	// 2 2 1
		 // ...
		 // sixth order moment
		 m(26, alpha) = e[alpha](0) * e[alpha](0) * e[alpha](1) * e[alpha](1)
		 * e[alpha](2) * e[alpha](2);	// 2 2 2
		 // revise !*/
	}
	//cout << "D3Q27:" << endl;
	//m.print(cout);
	return m;
}

}/* namespace natrium */

