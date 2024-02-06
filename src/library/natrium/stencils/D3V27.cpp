/*
 * D3V27.cpp

 *
 *  Created on: Oct 29, 2020
 *      Author: dwilde3m
 */

#include "D3V27.h"
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
D3V27::D3V27(double scaling) :
		Stencil(D, Q, makeDirections(scaling), makeWeights(), Stencil_D3V27,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
				scaling * pow(3, -0.5)), m_speedOfSoundSquare(
				scaling * scaling / 3.), m_scaling(scaling) {
} //constructor

/// destructor
D3V27::~D3V27() {
} /// destructor

// make weights
vector<double> D3V27::makeWeights() {
	const double w0 = 0.31247897198654906;//(720+8*sqrt(15))/2205;
	const double w1 = 0.029035130153906134;//(270-46*sqrt(15))/15435;
	const double w2 = 0.0005195469396656799;//(162+41*sqrt(15))/6174;
	const double w3 = 0.06338446047675325;//(783-202*sqrt(15))/24696;
    vector<double> result;
	result += w0,
			w1,
			w1,
			w1,
			w1,
			w1,
			w1,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w2,
			w3,
			w3,
			w3,
			w3,
			w3,
			w3,
			w3,
			w3;
	return result;
} /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
vector<numeric_vector> D3V27::makeDirections(double scaling) {
	const double r = 2.358709038202103;//sqrt((15+sqrt(15))/2.0);
    const double s = 3.142130383387586; //sqrt(6-sqrt(15));
    const double t = 1.1198362860638005;//sqrt(9+sqrt(15));


    const double directionsArray[][3] =
	{
			{ 0.0, 0.0, 0.0 },
			{ r, 0.0, 0.0 },
			{ -r, 0.0, 0.0 },
			{ 0.0, r, 0.0 },
			{ 0.0, -r, 0.0 },
			{ 0.0, 0.0, r },
			{ 0.0, 0.0, -r },
			{ 0.0, s,	s },
			{ 0.0, -s, -s },
			{ 0.0, s,	-s },
			{ 0.0, -s, s },
			{ s, 0.0,	s },
			{ -s, 0.0, -s },
			{ s, 0.0,	-s },
			{ -s, 0.0, s },
			{ s, s,	0.0 },
			{ -s, -s, 0.0 },
			{ s, -s, 0.0 },
			{ -s, s, 0.0 },
			{ t, t, t },
			{ -t, -t, -t },
			{t, t, -t },
			{ -t, -t, t },
			{ t, -t, t },
			{-t, t, -t },
			{ t, -t, -t },
			{ -t, t, t } };
	vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(D);
        direction(0) = scaling * directionsArray[i][0] / sqrt(3.);
        direction(1) = scaling * directionsArray[i][1] / sqrt(3.);
        direction(2) = scaling * directionsArray[i][2] / sqrt(3.);
		result += direction;
	}
	return result;
}
/// make directions

numeric_matrix D3V27::makeMomentBasis(vector<numeric_vector> e) {
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

}/* namespace natrium */

