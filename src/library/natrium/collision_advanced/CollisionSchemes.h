/*
 * CollisionSchemes.h
 *
 *  Created on: 19.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_
#define LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_

#include "Equilibria.h"
#include "SolverConfiguration.h"

namespace natrium {
template<int T_D, int T_Q, template<int T_D, int T_Q> class T_equilibrium>
class BGKCollision {
public:

	struct uniqueData {
		CollisionParameters<T_D, T_Q>& params;
		uniqueData(CollisionParameters<T_D, T_Q>& parameters) :
				params(parameters) {
			// empty for BGK collision
		}
	};

	void relax(double fLocal[], CollisionParameters<T_D, T_Q>& params,
			uniqueData& data) {
		double feq[T_Q];
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(feq, params);

		//Relax every direction towards the equilibrium
		for (int p = 0; p < T_Q; ++p) {
			fLocal[p] -= 1. / params.tau * (fLocal[p] - feq[p]);
		}
	}
};

template<int T_D, int T_Q, template<int T_D, int T_Q> class T_equilibrium>
class Regularized {
public:

	struct uniqueData {
		std::array<std::array<std::array<double, T_D>, T_D>, T_Q> Q = { { { } } };
		std::array<std::array<double, T_Q>, T_Q> pi = { { } };
		std::array<std::array<double, T_Q>, T_Q> pieq = { { } };

		CollisionParameters<T_D, T_Q>& params;

		uniqueData(CollisionParameters<T_D, T_Q>& parameters) :
				params(parameters) {
			initializeQ();
		}

		void initializeQ() {
			for (int a = 0; a < T_Q; a++) {
				for (int b = 0; b < T_D; b++) {
					for (int c = 0; c < T_D; c++) {
						Q[a][b][c] = params.e[a][b] * params.e[a][c];
						if (b == c) {
							Q[a][b][c] -= params.cs2;
						}

					}
				}
			}
		}

	};

	void relax(double fLocal[], CollisionParameters<T_D, T_Q>& params,
			uniqueData& data) {
		double feq[T_Q];
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(feq, params);

		for (int m = 0; m < T_D; m++) {
			for (int n = 0; n < T_D; n++) {
				data.pi[m][n] = 0.0;
				data.pieq[m][n] = 0.0;
			}
		}

		for (int j = 0; j < T_Q; j++) {
			for (int m = 0; m < T_D; m++) {
				for (int n = 0; n < T_D; n++) {
					data.pi[m][n] += fLocal[j] * params.e[j][m]
							* params.e[j][n];
					data.pieq[m][n] += feq[j] * params.e[j][m] * params.e[j][n];
				}
			}
		}



		for (int m = 0; m < T_D; m++) {
			for (int n = 0; n < T_D; n++) {
				data.pi[m][n] -= data.pieq[m][n];
			}
		}

		std::array<double, T_Q> fi1 = { 0.0 };

		for (int a = 0; a < T_Q; a++) {
			for (int b = 0; b < T_D; b++) {
				for (int c = 0; c < T_D; c++) {

					fi1[a] += params.weight[a] / (2 * params.cs2 * params.cs2)
							* data.Q[a][b][c] * data.pi[b][c];

				}
			}

		}

		//Relax every direction towards the equilibrium
		for (int i = 0; i < T_Q; ++i) {
			fLocal[i] = feq[i] + (1. - 1. / params.tau) * fi1[i];

		}

	} //relax
} //class regularized
;


template<int T_D, int T_Q, template<int T_D, int T_Q> class T_equilibrium>
class MRT {
public:

	struct uniqueData {
		const std::array<std::array<double, T_Q>, T_Q> M = { { } };
		const std::array<std::array<double, T_Q>, T_Q> T = { { } };
		const std::array<double, T_Q> omega = { { } };
		std::array<double, T_Q> m = { { } };
		std::array<double, T_Q> meq = { { } };
		CollisionParameters<T_D, T_Q>& params;

		uniqueData(CollisionParameters<T_D, T_Q>& parameters) :
				params(parameters) {
		}
	}; /* UniqueData */


	void relax(double fLocal[], CollisionParameters<T_D, T_Q>& params,
			uniqueData& data) {
		double feq[T_Q];
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(feq, params);

		// calculate moments
		AuxiliaryMRTFunctions::matrix_vector_product(data.M, fLocal, data.m);
		AuxiliaryMRTFunctions::matrix_vector_product(data.M, feq, data.meq);

		// relax moments
		for (size_t i = 0; i < T_Q; i++) {
			data.m[i] = data.m[i] - data.omega[i] * (data.m[i] - data.meq[i]);
		}
		// transform back
		AuxiliaryMRTFunctions::matrix_vector_product(data.T, data.m, fLocal);

	} //relax
} //class MRT
;
} // namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_ */
