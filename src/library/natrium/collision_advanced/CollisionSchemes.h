/*
 * CollisionSchemes.h
 *
 *  Created on: 19.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_
#define LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_

#include "Equilibria.h"
#include "../solver/SolverConfiguration.h"
#include "AuxiliaryMRTFunctions.h"

namespace natrium {
template<int T_D, int T_Q, template<int, int> class T_equilibrium>
class BGKCollision {
public:

	struct SpecificCollisionData {
		GeneralCollisionData<T_D, T_Q>& genData;
		SpecificCollisionData(GeneralCollisionData<T_D, T_Q>& parameters) :
				genData(parameters) {
			// empty for BGK collision
		}
	};

	void relax(std::array<double, T_Q>& fLocal,
			GeneralCollisionData<T_D, T_Q>& genData,
			SpecificCollisionData& specData) {
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(genData.feq, genData);

		//Relax every direction towards the equilibrium
		for (int p = 0; p < T_Q; ++p) {
			fLocal[p] -= 1. / genData.tau * (fLocal[p] - genData.feq[p]);
		}
	}

		void relaxWithG(std::array<double, T_Q>& fLocal, std::array<double, T_Q>& gLocal,
			GeneralCollisionData<T_D, T_Q>& genData,
			SpecificCollisionData& specData) {
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(genData.feq, genData);
		calculateGeqFromFeq<T_D,T_Q>(genData.feq,genData.geq,genData);

		//Distribution functions for variable Prandtl number (cf. Frapolli 2019)
		std::array<double, T_Q> fStar = genData.feq;
        std::array<double, T_Q> gStar = genData.geq;

        bool variablePrandtlNumber = true;

        if (variablePrandtlNumber == true) {

            double heatFluxTensorF[T_D][T_D][T_D] = {{{0.0}}};
            double heatFluxTensorFEq[T_D][T_D][T_D] = {{{0.0}}};
            double heatFluxTensorG[T_D][T_D][T_D] = {{{0.0}}};
            double heatFluxTensorGEq[T_D][T_D][T_D] = {{{0.0}}};

            calculateCenteredHeatFluxTensor<T_D,T_Q>(fLocal, heatFluxTensorF, genData);
            calculateCenteredHeatFluxTensor<T_D,T_Q>(genData.feq, heatFluxTensorFEq, genData);
            calculateCenteredHeatFluxTensor<T_D,T_Q>(gLocal, heatFluxTensorG, genData);
            calculateCenteredHeatFluxTensor<T_D,T_Q>(genData.geq,heatFluxTensorGEq,genData);

            calculateFStar<T_D, T_Q>(fStar, heatFluxTensorF, heatFluxTensorFEq, genData);
            calculateFStar<T_D, T_Q>(gStar, heatFluxTensorG, heatFluxTensorGEq, genData);
        }

        double sutherland_tau = genData.tau; // (genData.tau-0.5)*sqrt(genData.temperature)+0.5;
        double energy_tau=genData.tau;

        double prandtl = 0.7;
        double prandtl_tau = (genData.tau - 0.5)/prandtl + 0.5;

        //if(genData.maskShockSensor>0.5)
        //{
            //sutherland_tau = (1.0+sqrt(genData.maskShockSensor)*10.0)*sutherland_tau;
            //energy_tau=sutherland_tau;
        //}

		//Relax every direction towards the equilibrium
		for (int p = 0; p < T_Q; ++p) {
            fLocal[p] -= 1. / sutherland_tau * (fLocal[p] - genData.feq[p]) + (1./ sutherland_tau-1./prandtl_tau)*(fStar[p]-genData.feq[p]);
            gLocal[p] -= 1. / energy_tau * (gLocal[p] - genData.geq[p]) + (1./ sutherland_tau-1./prandtl_tau)*(gStar[p]-genData.feq[p]);;
		}
	}
};


template<int T_D, int T_Q, template<int, int > class T_equilibrium>
class Regularized {
public:

	struct SpecificCollisionData {
		std::array<std::array<std::array<double, T_D>, T_D>, T_Q> Q =
				{ { { } } };
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

	void relax(std::array<double, T_Q>& fLocal,
			GeneralCollisionData<T_D, T_Q>& genData,
			SpecificCollisionData& specData) {
		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(genData.feq, genData);

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
					specData.pieq[m][n] += genData.feq[j] * genData.e[j][m]
							* genData.e[j][n];
				}
			}
		}

		for (int m = 0; m < T_D; m++) {
			for (int n = 0; n < T_D; n++) {
				specData.pi[m][n] -= specData.pieq[m][n];
			}
		}

		std::array<double, T_Q> fi1 = {{ 0.0 }};

		for (int a = 0; a < T_Q; a++) {
			for (int b = 0; b < T_D; b++) {
				for (int c = 0; c < T_D; c++) {

					fi1[a] += genData.weight[a]
							/ (2 * genData.cs2 * genData.cs2)
							* specData.Q[a][b][c] * specData.pi[b][c];

				}
			}

		}

		//Relax every direction towards the equilibrium
		for (int i = 0; i < T_Q; ++i) {
			fLocal[i] = genData.feq[i] + (1. - 1. / genData.tau) * fi1[i];

		}

	} //relax
} //class regularized
;

template<int T_D, int T_Q, template<int, int> class T_equilibrium>
class MultipleRelaxationTime {
public:

	struct SpecificCollisionData {
		GeneralCollisionData<T_D, T_Q>& genData;
		const std::array<std::array<double, T_Q>, T_Q> M;
		const std::array<std::array<double, T_Q>, T_Q> T;
		const std::array<double, T_Q> omega;
		std::array<double, T_Q> m = { { } };
		std::array<double, T_Q> meq = { { } };

		SpecificCollisionData(GeneralCollisionData<T_D, T_Q>& genData) :
				genData(genData), M(
						AuxiliaryMRTFunctions::make_M<T_Q>(
								genData.configuration.getMRTBasis())), T(
						AuxiliaryMRTFunctions::make_T<T_Q>(
								genData.configuration.getMRTBasis())), omega(
						AuxiliaryMRTFunctions::make_diag<T_Q>(genData.tau,
								genData.configuration.getMRTBasis(),
								genData.configuration.getMRTRelaxationTimes())) {
		}
	}; /* SpecificCollisionData */

	void relax(std::array<double, T_Q>& fLocal,
			GeneralCollisionData<T_D, T_Q>& genData,
			SpecificCollisionData& specData) {

		//Initialize the corresponding Equilibrium Distribution Function
		T_equilibrium<T_D, T_Q> eq;

		//Calculate the equilibrium and write the result to feq
		eq.calc(genData.feq, genData);

		// calculate moments
		AuxiliaryMRTFunctions::matrix_vector_product(specData.M, fLocal,
				specData.m);
		AuxiliaryMRTFunctions::matrix_vector_product(specData.M, genData.feq,
				specData.meq);

		// relax moments
		for (size_t i = 0; i < T_Q; i++) {
			specData.m[i] = specData.m[i]
					- specData.omega[i] * (specData.m[i] - specData.meq[i]);
		}

		// transform back
		AuxiliaryMRTFunctions::matrix_vector_product(specData.T, specData.m,
				fLocal);

	} //relax
} //class MRT
;

} // namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_COLLISIONSCHEMES_H_ */
