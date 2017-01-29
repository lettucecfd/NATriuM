/*
 * CollisionSchemes.h
 *
 *  Created on: 19.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_
#define LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_

#include "Equilibria.h"
#include "AuxiliaryCollisionFunctions.h"

namespace natrium{
template<int T_D, int T_Q, template <int T_D, int T_Q> class T_equilibrium>
class BGKCollision {
public:

	void relax(double fLocal[], CollisionParameters<T_D,T_Q>& params) {
		double feq[T_Q];
		T_equilibrium<T_D,T_Q> eq;

		eq.calc(feq, params);
		cout << calculateDensity<9>(feq) << " " << calculateDensity<9>(fLocal) <<  endl;

		for (int p = 0; p < T_Q; ++p) {
		fLocal[p]= fLocal[p] - 1./params.tau*(fLocal[p]-feq[p]);
	}
}
};
}// namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_ */
