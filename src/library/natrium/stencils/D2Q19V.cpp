/**
 * @file D2Q19V.cpp
 * @short 
 * @date 02.06.2020
 * @author dwilde3m, University of Siegen, Germany
 */

#include "D2Q19V.h"

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
D2Q19V::D2Q19V(double scaling) :
		Stencil(D, Q, makeDirections(scaling), makeWeights(), Stencil_D2Q19V,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
                scaling/sqrt(3)), m_speedOfSoundSquare(
                scaling * scaling/(3)), m_scaling(scaling) {
} //constructor

/// destructor
D2Q19V::~D2Q19V() {
} /// destructor



// make weights
vector<double> D2Q19V::makeWeights() {
//    double r = (sqrt(5.)-sqrt(2.))/sqrt(3.);
//    double w_0 = (-3-3*r*r*r*r+54*r*r)/(75*r*r);
//    double w_m = (9*r*r*r*r-6-27*r*r)/(300*r*r*(r*r-1));
//    double w_n = (9-6*r*r*r*r-27*r*r)/(300*(1-r*r));
//    double w_0n = w_0*w_n;
//    double w_0m = w_0*w_m;
//    double w_mm = w_m*w_m;
//    double w_mn = w_m*w_n;
//    double w_nn = w_n*w_n;
	vector<double> result { 0.3168437267921905 , 0.1024247123210936 ,
                           0.1024247123210936 ,
                           0.1024247123210936 ,
                           0.1024247123210936 ,
                           0.002405335328939458 ,
                           0.002405335328939458 ,
                           0.002405335328939458 ,
                           0.002405335328939458 ,
                           0.00953510698543825 ,
                           0.00953510698543825 ,
                           0.00953510698543825 ,
                           0.00953510698543825 ,
                           0.006865104210104631 ,
                           0.006865104210104631 ,
                           0.10558878375062891 ,
                           0.10558878375062891 ,
                           0.0003939393722285871 ,
                           0.0003939393722285871
                          };
	//result += w_0*w_0, w_0m,w_0m, w_0m,w_0m, w_mm, w_mm, w_mm, w_mm, w_0n, w_0n, w_0n, w_0n, w_nn, w_nn, w_nn, w_nn, w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn ;
	return result;
} /// make weights

/// make directions
vector<numeric_vector> D2Q19V::makeDirections(double scaling) {

        vector<numeric_vector> result;

        const double directionsArray[19][2] {{ 0.0 , 0.0 },{ 1.367469636752619 , 0.775196278121181 },
            { 1.367469636752619 , -0.775196278121181 },
            { -1.367469636752619 , 0.775196278121181 },
            { -1.367469636752619 , -0.775196278121181 },
            { 2.6987507639352253 , 1.8663975507141328 },
            { 2.6987507639352253 , -1.8663975507141328 },
            { -2.6987507639352253 , 1.8663975507141328 },
            { -2.6987507639352253 , -1.8663975507141328 },
            { 1.105629214668943 , 2.5175897644357486 },
            { 1.105629214668943 , -2.5175897644357486 },
            { -1.105629214668943 , 2.5175897644357486 },
            { -1.105629214668943 , -2.5175897644357486 },
            { 2.9213306655318734 , 0.0 },
            { -2.9213306655318734 , 0.0 },
            { 0.0 , 1.4869982213169028 },
            { 0.0 , -1.4869982213169028 },
            { 0.0 , 3.8358342053914734 },
            { 0.0 , -3.8358342053914734 }};

	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(2);
		direction(0) = scaling * directionsArray[i][0] / sqrt(3);
		direction(1) = scaling * directionsArray[i][1] / sqrt(3);
		result += direction;
	}
	return result;
} /// make directions

numeric_matrix D2Q19V::makeMomentBasis(vector<numeric_vector> e) {
    (void)e;
    numeric_matrix m(Q);
    for (size_t i = 0;i<Q;i++){
        for (size_t j = 0;j<Q;j++){
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

