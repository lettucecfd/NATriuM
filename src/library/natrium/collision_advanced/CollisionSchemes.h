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
		std::array<double, T_Q> fStar = {0.0};//genData.feq;
        std::array<double, T_Q> gStar = {0.0};//genData.geq;

        std::array<double, T_Q> fNeq = {0.0};
        std::array<double, T_Q> gNeq = {0.0};
#pragma omp parallel for
        for (int p = 0; p < T_Q; ++p) {
            fNeq[p] = fLocal[p] - genData.feq[p];
            gNeq[p] = gLocal[p] - genData.geq[p];
        }

        const bool isPrandtlNumberSet = genData.configuration.isPrandtlNumberSet();

        if (isPrandtlNumberSet == true) {

            // 3 staged non-equilibrium heat flux tensors
            std::array<std::array<std::array<double, T_D>, T_D>, T_D> heatFluxTensorFNEq = {{{0.0}}};
            std::array<std::array<std::array<double, T_D>, T_D>, T_D> heatFluxTensorGNeq = {{{0.0}}};

            calculateCenteredHeatFluxTensor<T_D,T_Q>(fLocal, heatFluxTensorFNEq, genData);
            //calculateCenteredHeatFluxTensor<T_D,T_Q>(genData.feq, heatFluxTensorFEq, genData);
            calculateCenteredHeatFluxTensor<T_D,T_Q>(gLocal, heatFluxTensorGNeq, genData);
            //calculateCenteredHeatFluxTensor<T_D,T_Q>(genData.geq,heatFluxTensorGEq,genData);

            calculateFStar<T_D, T_Q>(fStar, heatFluxTensorFNEq, genData);
            calculateFStar<T_D, T_Q>(gStar, heatFluxTensorGNeq, genData);
        }
        double sutherland_factor = 1.402*pow(genData.temperature,1.5) / ( genData.temperature + 0.40417);
        double visc_tau = (genData.tau-0.5)*sutherland_factor/(genData.temperature*genData.density)+0.5;


        const double knudsen_estimate = calculateKnudsenNumberEstimate<T_D, T_Q>(fLocal, genData.feq, genData.weight);
        double tau_factor = 1.0;
            if(knudsen_estimate >= 0.01)
                tau_factor = 1.05;
            if(knudsen_estimate >= 0.05)
                tau_factor = 1.35;
            if(knudsen_estimate >= 0.1)
                tau_factor = 1/visc_tau;
        visc_tau *=tau_factor;
        double ener_tau = visc_tau;
        double prandtl = genData.configuration.getPrandtlNumber();
        double prandtl_tau = (visc_tau - 0.5) / prandtl + 0.5;

        genData.maskShockSensor = knudsen_estimate;

        //if(genData.maskShockSensor>0.5)
        //{
            //visc_tau = (1.0+sqrt(genData.maskShockSensor)*10.0)*visc_tau;
            //ener_tau=visc_tau;
        //}

		//Relax every direction towards the equilibrium
#pragma omp parallel for
            for (int p = 0; p < T_Q; ++p) {
            fLocal[p] -= 1. / visc_tau * (fNeq[p]) + (1. / visc_tau - 1. / prandtl_tau) * fStar[p]; // -genData.feq[p]);
            gLocal[p] -= 1. / ener_tau * (gNeq[p]) + (1. / visc_tau - 1. / prandtl_tau) * gStar[p];
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
