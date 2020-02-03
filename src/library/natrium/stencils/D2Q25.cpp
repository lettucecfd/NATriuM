/**
 * @file D2Q25.cpp
 * @short 
 * @date 09.10.2018
 * @author dwilde3m, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "D2Q25.h"

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
const size_t D2Q25::D = 2;
/// Q
const size_t D2Q25::Q = 25;

/// constructor
D2Q25::D2Q25(double scaling) :
		Stencil(2, 25, makeDirections(scaling), makeWeights(), Stencil_D2Q25,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
                scaling/sqrt(3.)), m_speedOfSoundSquare(
                scaling * scaling/(3.)), m_scaling(scaling) {
} //constructor

/// destructor
D2Q25::~D2Q25() {
} /// destructor



// make weights
vector<double> D2Q25::makeWeights() {
	double w_0 = 4./45*(4.+sqrt(10.));
	double w_m = 3./80*(8.-sqrt(10.));
	double w_n = 1./720*(16.-5.*sqrt(10.));


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
vector<numeric_vector> D2Q25::makeDirections(double scaling) {
    const double c_m = scaling * sqrt(1. - sqrt(2./5)) /sqrt(3.);
    const double c_n = 3.0*scaling* sqrt(1. - sqrt(2./5)) /sqrt(3.);
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

numeric_matrix D2Q25::makeMomentBasis(vector<numeric_vector> e) {
	numeric_matrix m(Q);
	for (int i = 0;i<Q;i++){
        for (int j = 0;j<Q;j++){
            if(i==j)
            {
               m(i,j) =1.0;
            } else
            {
                m(i,j) =0.0;
            }
        }

	}
	return m;
}

} /* namespace natrium */

