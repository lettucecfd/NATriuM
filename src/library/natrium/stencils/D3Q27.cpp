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
		Stencil(3, 27, makeDirections(scaling), makeWeights(), Stencil_D3Q27), m_speedOfSound(
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
			2. / 27., 2. / 27., 2. / 27., 2. / 27., 2. / 27., 2. / 27.,
			1. / 54., 1. / 54., 1. / 54., 1. / 54., 1. / 54., 1. / 54.,
			1. / 54., 1. / 54., 1. / 54., 1. / 54., 1. / 54., 1. / 54.,
			1. / 216., 1. / 216., 1. / 216., 1. / 216., 1.
			/ 216., 1. / 216., 1. / 216., 1. / 216.;
	return result;
} /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
vector<numeric_vector> D3Q27::makeDirections(double scaling) {
	const double directionsArray[][3] = { { 0.0, 0.0, 0.0 },
			{ scaling, 0.0, 0.0 }, { -scaling, 0.0, 0.0 },
			{ 0.0, scaling, 0.0 }, { 0.0, -scaling, 0.0 },
			{ 0.0, 0.0, scaling }, { 0.0, 0.0, -scaling }, { 0.0, scaling,
					scaling }, { 0.0, -scaling - scaling, }, { 0.0, scaling,
					-scaling }, { 0.0, -scaling, scaling }, { scaling, 0.0,
					scaling }, { -scaling, 0.0, -scaling }, { scaling, 0.0,
					-scaling }, { -scaling, 0.0, scaling }, { scaling, scaling,
					0.0 }, { -scaling, -scaling, 0.0 },
			{ scaling, -scaling, 0.0 }, { -scaling, scaling, 0.0}, { scaling,
					scaling, scaling }, { -scaling, -scaling, -scaling }, {
					scaling, scaling, -scaling },
			{ -scaling, -scaling, scaling }, { scaling, -scaling, scaling }, {
					-scaling, scaling, -scaling },
			{ scaling, -scaling, -scaling }, { -scaling, scaling, scaling } };
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

}/* namespace natrium */

