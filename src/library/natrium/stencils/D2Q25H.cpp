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
const size_t D2Q25H::Q = 19;

/// constructor
D2Q25H::D2Q25H(double scaling) :
		Stencil(2, 19, makeDirections(scaling), makeWeights(), Stencil_D2Q25H,
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

