/*
 * Equilibria.h
 *
 *  Created on: 19.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_
#define LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_
//#include "CollisionOperator.h"
#include "AuxiliaryCollisionFunctions.h"
#include "Equilibria.h"
#include "CollisionSchemes.h"
namespace natrium{
template <int T_D,int T_Q>
class BGKEquilibrium
{

public:
	//BGKEquilibrium();
	void calc(double feq[], const CollisionParameters<T_D,T_Q> & params);
};


template<>
 inline void BGKEquilibrium<2,9>::calc(double feq[], const CollisionParameters<2,9> & params) {
	double scalar_product;
	double weighting;
	double uSquareTerm;
	double mixedTerm;
	double prefactor = params.scaling / params.cs2;
	// calculate equilibrium distribution
	scalar_product = params.velocity[0] * params.velocity[0]
			+ params.velocity[1] * params.velocity[1];
	uSquareTerm = -scalar_product / (2 * params.cs2);
	// direction 0
	weighting = 4. / 9. * params.density;
	feq[0] = weighting * (1 + uSquareTerm);
	// directions 1-4
	weighting = 1. / 9. * params.density;
	mixedTerm = prefactor * (params.velocity[0]);
	feq[1] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
	feq[3] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
	mixedTerm = prefactor * (params.velocity[1]);
	feq[2] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
	feq[4] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
	// directions 5-8
	weighting = 1. / 36. * params.density;
	mixedTerm = prefactor * (params.velocity[0] + params.velocity[1]);
	feq[5] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
	feq[7] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
	mixedTerm = prefactor * (-params.velocity[0] + params.velocity[1]);
	feq[6] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
	feq[8] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
}
}

#endif /* LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_ */
