/*
 * CollisionSchemes.h
 *
 *  Created on: 19.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_
#define LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_

#include "Equilibria.h"

namespace natrium {
template<int T_D, int T_Q, template<int T_D, int T_Q> class T_equilibrium>
class BGKCollision {
public:

	struct SpecificCollisionData {
		GeneralCollisionData<T_D, T_Q>& genData;
		SpecificCollisionData(GeneralCollisionData<T_D, T_Q>& parameters) :
				genData(parameters) {
			// empty for BGK collision
		}
	};

	void relax(std::array<double,T_Q>& fLocal, GeneralCollisionData<T_D, T_Q>& genData,
			SpecificCollisionData& specData) {
		double feq[T_Q];
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(feq, genData);

		//Relax every direction towards the equilibrium
		for (int p = 0; p < T_Q; ++p) {
			fLocal[p] -= 1. / genData.tau * (fLocal[p] - feq[p]);
		}
	}
};

template<int T_D, int T_Q, template<int T_D, int T_Q> class T_equilibrium>
class Regularized {
public:

	struct SpecificCollisionData {
		std::array<std::array<std::array<double, T_D>, T_D>, T_Q> Q = { { { } } };
		std::array<std::array<double, T_Q>, T_Q> pi = { { } };
		std::array<std::array<double, T_Q>, T_Q> pieq = { { } };

		GeneralCollisionData<T_D, T_Q>& genData;

		SpecificCollisionData(GeneralCollisionData<T_D, T_Q>& parameters) :
				genData(parameters) {
			initializeQ();
		}

		void initializeQ() {
			for (int a = 0; a < T_Q; a++) {
				for (int b = 0; b < T_D; b++) {
					for (int c = 0; c < T_D; c++) {
						Q[a][b][c] = genData.e[a][b] * genData.e[a][c];
						if (b == c) {
							Q[a][b][c] -= genData.cs2;
						}

					}
				}
			}
		}

	};

	void relax(std::array<double,T_Q>& fLocal, GeneralCollisionData<T_D, T_Q>& genData,
			SpecificCollisionData& specData) {
		double feq[T_Q];
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(feq, genData);

		for (int m = 0; m < T_D; m++) {
			for (int n = 0; n < T_D; n++) {
				specData.pi[m][n] = 0.0;
				specData.pieq[m][n] = 0.0;
			}
		}

		for (int j = 0; j < T_Q; j++) {
			for (int m = 0; m < T_D; m++) {
				for (int n = 0; n < T_D; n++) {
					specData.pi[m][n] += fLocal[j] * genData.e[j][m]
							* genData.e[j][n];
					specData.pieq[m][n] += feq[j] * genData.e[j][m] * genData.e[j][n];
				}
			}
		}



		for (int m = 0; m < T_D; m++) {
			for (int n = 0; n < T_D; n++) {
				specData.pi[m][n] -= specData.pieq[m][n];
			}
		}

		std::array<double, T_Q> fi1 = { 0.0 };

		for (int a = 0; a < T_Q; a++) {
			for (int b = 0; b < T_D; b++) {
				for (int c = 0; c < T_D; c++) {

					fi1[a] += genData.weight[a] / (2 * genData.cs2 * genData.cs2)
							* specData.Q[a][b][c] * specData.pi[b][c];

				}
			}

		}

		//Relax every direction towards the equilibrium
		for (int i = 0; i < T_Q; ++i) {
			fLocal[i] = feq[i] + (1. - 1. / genData.tau) * fi1[i];

		}

	} //relax
} //class regularized
;
} // namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_ */
