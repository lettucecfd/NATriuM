/*
 * RD3Q27.cpp

 *
 *  Created on: Nov 20, 2015
 *      Author: kraemer
 */

#include "RD3Q27.h"
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
RD3Q27::RD3Q27(double scaling) :
		Stencil(3, 27, makeDirections(scaling), makeWeights(), Stencil_RD3Q27,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
				scaling * pow(5, -0.5)), m_speedOfSoundSquare(
				scaling * scaling / 5.), m_scaling(scaling) {
} //constructor

/// destructor
RD3Q27::~RD3Q27() {
} /// destructor

// make weights
vector<double> RD3Q27::makeWeights() {
	vector<double> result;
	result += 1. / 3.,
			1. / 30.,
			1. / 30.,
			1. / 30.,
			1. / 30.,
			1. / 30.,
			1. / 30.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			1. / 300.,
			4. / 75.,
			4. / 75.,
			4. / 75.,
			4. / 75.,
			4. / 75.,
			4. / 75.,
			4. / 75.,
			4. / 75.;
	return result;
} /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
vector<numeric_vector> RD3Q27::makeDirections(double scaling) {
	double halfScaling = scaling / 2.;
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
			{ halfScaling, halfScaling, halfScaling },
			{ -halfScaling, -halfScaling, -halfScaling },
			{halfScaling, halfScaling, -halfScaling },
			{ -halfScaling, -halfScaling, halfScaling },
			{ halfScaling, -halfScaling, halfScaling },
			{-halfScaling, halfScaling, -halfScaling },
			{ halfScaling, -halfScaling, -halfScaling },
			{ -halfScaling, halfScaling, halfScaling } };
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

numeric_matrix RD3Q27::makeMomentBasis(vector<numeric_vector> e) {
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

