/**
 * @file D2Q9Model.cpp
 * @short 
 * @date 30.08.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "D2Q9Model.h"

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
const size_t D2Q9Model::D = 2;
/// Q
const size_t D2Q9Model::Q = 9;

/// constructor
D2Q9Model::D2Q9Model(double scaling):
		BoltzmannModel(2, 9, makeDirections(scaling), makeWeights(), Stencil_D2Q9),
		m_speedOfSound(scaling*pow(3, -0.5)),
		m_speedOfSoundSquare(scaling*scaling/3.),
		m_scaling(scaling){
} //constructor


/// destructor
D2Q9Model::~D2Q9Model() {
} /// destructor


// make weights
vector<double> D2Q9Model::makeWeights()  {
	vector<double> result;
	result += 4./9., 1./9., 1./9., 1./9., 1./9.,
			1./36., 1./36., 1./36., 1./36.;
	return result;
}/// make weights


/// make directions
vector<numeric_vector> D2Q9Model::makeDirections(double scaling) {
	const double directionsArray[][2] = { { 0.0, 0.0 }, { scaling, 0.0 }, { 0.0,
			scaling }, { -scaling, 0.0 }, { 0.0, -scaling }, { scaling, scaling },
			{ -scaling, scaling }, { -scaling, -scaling }, { scaling, -scaling } };
	vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(2);
		direction(0) = directionsArray[i][0];
		direction(1) = directionsArray[i][1];
		result += direction;
	}
	return result;
}/// make directions



} /* namespace natrium */
