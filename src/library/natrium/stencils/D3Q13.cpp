/*
 * D3Q13.cpp

 *
 *  Created on: Oct 20, 2017
 *      Author: dwilde
 */

#include "D3Q13.h"
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
    const size_t D3Q13::D = 3;
/// Q
    const size_t D3Q13::Q = 13;

/// constructor
    D3Q13::D3Q13(double scaling) :
            Stencil(3, 13, makeDirections(scaling), makeWeights(), Stencil_D3Q13,
                    makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
            scaling * pow(3, -0.5)), m_speedOfSoundSquare(
            scaling * scaling / 3.), m_scaling(scaling) {
        if (scaling > 100) {
            LOG(WARNING)
                    << "The D3Q13 stencil is used with scaling > 100. "
                       "This may lead to significant round-off errors."
                       "(See UnitTest for Moment matrix of D3Q13 model for details.)"
                    << endl;
        }
    } //constructor

/// destructor
    D3Q13::~D3Q13() {
    } /// destructor

// make weights
    vector<double> D3Q13::makeWeights() {
        vector<double> result;
        result += 0.4, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05;
        return result;
    } /// make weights

/// make directions
    vector <numeric_vector> D3Q13::makeDirections(double scaling) {
        // Scale to the standard speed of sound
        scaling *= 1./sqrt(3);


        const double phi_plus = scaling * sqrt((1. + sqrt(5)) / 2.);
        const double phi_minu = scaling * sqrt((1. - sqrt(5)) / 2.);
        const double directionsArray[][3] = {{0.0,       0.0,       0.0},
                                             {0.0,       phi_minu,  phi_plus},
                                             {0.0,       -phi_minu, -phi_plus},
                                             {0.0,       -phi_minu, phi_plus},
                                             {0.0,       phi_minu,  -phi_plus},
                                             {phi_minu,  phi_plus,  0.0},
                                             {-phi_minu, -phi_plus, 0.0},
                                             {-phi_minu, phi_plus,  0.0},
                                             {phi_minu,  -phi_plus, 0.0},
                                             {phi_plus,  0.0,       phi_minu},
                                             {-phi_plus, 0.0,       -phi_minu},
                                             {-phi_plus, 0.0,       phi_minu},
                                             {phi_plus,  0.0,       -phi_minu}};
        vector<numeric_vector> result;
        for (size_t i = 0; i < Q; i++) {
            numeric_vector direction(D);
            direction(0) = directionsArray[i][0];
            direction(1) = directionsArray[i][1];
            direction(2) = directionsArray[i][2];
            result += direction;
        }
        return result;
    } /// make directions

    numeric_matrix D3Q13::makeMomentBasis(vector <numeric_vector> e) {
        // This function is without any use in D3Q13
        numeric_matrix m(Q);
        return m;
    }

}/* namespace natrium */

