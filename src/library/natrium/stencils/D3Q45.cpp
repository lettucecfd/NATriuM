/**
 * @file D3Q45.cpp
 * @short 
 * @date 02.06.2020
 * @author dwilde3m, University of Siegen, Germany
 */

#include "D3Q45.h"

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
D3Q45::D3Q45(double scaling) :
		Stencil(D, Q, makeDirections(scaling), makeWeights(), Stencil_D3Q45,
				makeMomentBasis(makeDirections(scaling))), m_speedOfSound(
                scaling/sqrt(3)), m_speedOfSoundSquare(
                scaling * scaling/(3)), m_scaling(scaling) {
} //constructor

/// destructor
D3Q45::~D3Q45() {
} /// destructor



// make weights
vector<double> D3Q45::makeWeights() {

	vector<double> result {
            0.20740740740740618 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.05787037037037047 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.00462962962962958 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
            0.0004629629629629939 ,
    };
	//result += w_0*w_0, w_0m,w_0m, w_0m,w_0m, w_mm, w_mm, w_mm, w_mm, w_0n, w_0n, w_0n, w_0n, w_nn, w_nn, w_nn, w_nn, w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn,w_mn ;
	return result;
} /// make weights

/// make directions
vector<numeric_vector> D3Q45::makeDirections(double scaling) {

        vector<numeric_vector> result;
	for (size_t i = 0; i < Q; i++) {
		numeric_vector direction(3);

        const double m_directionsArray[45][3] = {
                {0.0 ,  0.0 ,  0.0},
                {0.06386083877343968 ,  -1.2239121278243665 ,  -1.2239121278243665},
                {-1.2239121278243665 ,  0.06386083877343968 ,  -1.2239121278243665},
                {-1.2239121278243665 ,  -1.2239121278243665 ,  0.06386083877343968},
                {1.5766994272507744 ,  -0.5069610024977665 ,  -0.5069610024977665},
                {-0.5069610024977665 ,  1.5766994272507744 ,  -0.5069610024977665},
                {-0.5069610024977665 ,  -0.5069610024977665 ,  1.5766994272507744},
                {0.5069610024977665 ,  0.5069610024977665 ,  -1.5766994272507744},
                {0.5069610024977665 ,  -1.5766994272507744 ,  0.5069610024977665},
                {-1.5766994272507744 ,  0.5069610024977665 ,  0.5069610024977665},
                {1.2239121278243665 ,  1.2239121278243665 ,  -0.06386083877343968},
                {1.2239121278243665 ,  -0.06386083877343968 ,  1.2239121278243665},
                {-0.06386083877343968 ,  1.2239121278243665 ,  1.2239121278243665},
                {2.9239876105912574 ,  0.4744978678080795 ,  0.4744978678080795},
                {0.4744978678080795 ,  2.9239876105912574 ,  0.4744978678080795},
                {0.4744978678080795 ,  0.4744978678080795 ,  2.9239876105912574},
                {1.7320508075688787 ,  1.7320508075688787 ,  1.7320508075688787},
                {2.403092127540177 ,  0.8892242114059369 ,  -1.5602655313772367},
                {2.403092127540177 ,  -1.5602655313772367 ,  0.8892242114059369},
                {1.5602655313772367 ,  -0.8892242114059369 ,  -2.403092127540177},
                {1.5602655313772367 ,  -2.403092127540177 ,  -0.8892242114059369},
                {0.8892242114059369 ,  2.403092127540177 ,  -1.5602655313772367},
                {0.8892242114059369 ,  -1.5602655313772367 ,  2.403092127540177},
                {-0.8892242114059369 ,  1.5602655313772367 ,  -2.403092127540177},
                {-0.8892242114059369 ,  -2.403092127540177 ,  1.5602655313772367},
                {-1.5602655313772367 ,  2.403092127540177 ,  0.8892242114059369},
                {-1.5602655313772367 ,  0.8892242114059369 ,  2.403092127540177},
                {-2.403092127540177 ,  1.5602655313772367 ,  -0.8892242114059369},
                {-2.403092127540177 ,  -0.8892242114059369 ,  1.5602655313772367},
                {-1.7320508075688787 ,  -1.7320508075688787 ,  -1.7320508075688787},
                {-0.4744978678080795 ,  -0.4744978678080795 ,  -2.9239876105912574},
                {-0.4744978678080795 ,  -2.9239876105912574 ,  -0.4744978678080795},
                {-2.9239876105912574 ,  -0.4744978678080795 ,  -0.4744978678080795},
                {2.7367507163016924 ,  2.7367507163016924 ,  -0.14279717659756475},
                {2.7367507163016924 ,  -0.14279717659756475 ,  2.7367507163016924},
                {-0.14279717659756475 ,  2.7367507163016924 ,  2.7367507163016924},
                {-3.5256070994177073 ,  1.1335992635264445 ,  1.1335992635264445},
                {1.1335992635264445 ,  -3.5256070994177073 ,  1.1335992635264445},
                {1.1335992635264445 ,  1.1335992635264445 ,  -3.5256070994177073},
                {-1.1335992635264445 ,  -1.1335992635264445 ,  3.5256070994177073},
                {-1.1335992635264445 ,  3.5256070994177073 ,  -1.1335992635264445},
                {3.5256070994177073 ,  -1.1335992635264445 ,  -1.1335992635264445},
                {0.14279717659756475 ,  -2.7367507163016924 ,  -2.7367507163016924},
                {-2.7367507163016924 ,  0.14279717659756475 ,  -2.7367507163016924},
                {-2.7367507163016924 ,  -2.7367507163016924 ,  0.14279717659756475}
        };


        direction(0) = scaling * m_directionsArray[i][0] / sqrt(3);
		direction(1) = scaling * m_directionsArray[i][1] / sqrt(3);
        direction(2) = scaling * m_directionsArray[i][2] / sqrt(3);
		result += direction;
	}
	return result;
} /// make directions

    numeric_matrix D3Q45::makeMomentBasis(vector<numeric_vector> e) {
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
