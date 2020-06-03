/**
 * @file D3Q77.cpp
 * @short 
 * @date 02.06.2020
 * @author dwilde3m, University of Siegen, Germany
 */

#include "D3Q77.h"

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
const size_t D3Q77::D = 3;
/// Q
const size_t D3Q77::Q = 77;

/// constructor
D3Q77::D3Q77(double scaling) :
		Stencil(3, 77, makeDirections(scaling), makeWeights(), Stencil_D3Q77,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
                scaling/sqrt(3)), m_speedOfSoundSquare(
                scaling * scaling/(3)), m_scaling(scaling) {
} //constructor

/// destructor
D3Q77::~D3Q77() {
} /// destructor



// make weights
vector<double> D3Q77::makeWeights() {

	vector<double> result {
            0.12148148148148114 ,
            0.0009194665015833565 ,
            0.0009194665015833565 ,
            0.0009194665015833565 ,
            0.0009194665015833565 ,
            0.0009194665015833565 ,
            0.0009194665015833565 ,
            0.08056201497989846 ,
            0.08056201497989846 ,
            0.08056201497989846 ,
            0.08056201497989846 ,
            0.08056201497989846 ,
            0.08056201497989846 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.016439238378214215 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            0.0025000000000000053 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            4.224310326716467e-05 ,
            0.01643923837821427 ,
            0.01643923837821427 ,
            0.01643923837821427 ,
            0.01643923837821427 ,
            0.01643923837821427 ,
            0.01643923837821427 ,
            0.01643923837821427 ,
            0.01643923837821427
    };
	//result += w_0*w_0, w_0m,w_0m, w_0m,w_0m, w_mm, w_mm, w_mm, w_mm, w_0n, w_0n, w_0n, w_0n, w_nn, w_nn, w_nn, w_nn, w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn ;
	return result;
} /// make weights

/// make directions
vector<numeric_vector> D3Q77::makeDirections(double scaling) {

        vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(3);
		direction(0) = scaling * m_directionsArray[i][0] / sqrt(3);
		direction(1) = scaling * m_directionsArray[i][1] / sqrt(3);
        direction(2) = scaling * m_directionsArray[i][2] / sqrt(3);
		result += direction;
	}
	return result;
} /// make directions

    numeric_matrix D3Q77::makeMomentBasis(vector<numeric_vector> e) {
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

