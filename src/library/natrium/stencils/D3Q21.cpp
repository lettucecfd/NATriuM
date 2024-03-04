/*
 * D3Q21.cpp

 *
 *  Created on: Oct 22, 2017
 *      Author: dwilde
 */

#include "D3Q21.h"
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
D3Q21::D3Q21(double scaling) :
		Stencil(D, Q, makeDirections(scaling), makeWeights(), Stencil_D3Q21,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
        scaling * pow(3, -0.5)), m_speedOfSoundSquare(
        scaling * scaling / 3.), m_scaling(scaling) {
} //constructor

/// destructor
D3Q21::~D3Q21() {
} /// destructor

// make weights
    vector<double> D3Q21::makeWeights() {
        vector<double> result;
        result += 2. / 5.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.,
                3. / 100.;
        return result;
    } /// make weights

/// make directions
/* Following the definition by Mohamad LBM book */
    vector <numeric_vector> D3Q21::makeDirections(double scaling) {
        //Scale by the nominal speed of sound of D3Q13 (0.4472) to obtain unity speed of sound
        scaling /= 0.774596669241483;
        // Then scale to the standard speed of sound
        scaling *= 1./sqrt(3);
        double phi = scaling * (1. + sqrt(5.)) / 2.;
        double inphi = scaling * 1. / ((1. + sqrt(5.)) / 2.);

        const double directionsArray[][3] =
                {
                        {0.0,      0.0,      0.0},
                        {scaling,  scaling,  scaling},
                        {-scaling, -scaling, -scaling},
                        {-scaling, scaling,  scaling},
                        {scaling,  -scaling, -scaling},
                        {scaling,  -scaling, scaling},
                        {-scaling, scaling,  -scaling},
                        {scaling,  scaling,  -scaling},
                        {-scaling, -scaling, scaling},
                        {0.0,      phi,      inphi},
                        {0.0,      -phi,     -inphi},
                        {0.0,      -phi,     inphi},
                        {0.0,      phi,      -inphi},
                        {inphi,    0.0,      phi},
                        {-inphi,   0.0,      -phi},
                        {-inphi,   0.0,      phi},
                        {inphi,    0.0,      -phi},
                        {phi,      inphi,    0.0},
                        {-phi,     -inphi,   0.0},
                        {-phi,     inphi,    0.0},
                        {phi,      -inphi,   0.0}};

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

numeric_matrix D3Q21::makeMomentBasis(vector<numeric_vector> e) {
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

