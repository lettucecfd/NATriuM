/*
 * AuxiliaryCollisionFunctions.h
 *
 *  Created on: 16.01.2017
 *      Author: natrium
 */
#ifndef AUXILIARYCOLLISIONFUNCTIONS_H_
#define AUXILIARYCOLLISIONFUNCTIONS_H_

#include "../collision/ExternalForceFunctions.h"
#include <array>
#include <vector>

namespace natrium {

template<int T_Q>
inline double calculateDensity(const double fLocal[]) {
	double density = 0;
	for (int p = 0; p < T_Q; ++p) {
		density += fLocal[p];
	}
	return density;
}

inline double calculateTauFromNu(double viscosity, double cs2,
		double timeStepSize) {
	double tau;
	tau = (viscosity) / (timeStepSize * cs2) + 0.5;
	return tau;
}

template<int T_Q>
inline void copyGlobalToLocalF(double fLocal[], const DistributionFunctions& f,
		const size_t i) {
	for (int p = 0; p < T_Q; ++p) {
		fLocal[p] = f.at(p)(i);
	}
}

template<int T_Q>
inline void copyLocalToGlobalF(double fLocal[], DistributionFunctions& f,
		const size_t i) {
	for (int p = 0; p < T_Q; ++p) {
		distributed_vector& f0 = f.at(p);
		f0(i) = fLocal[p];
	}
}
//Stores the needed parameters for the collision phase
template<int T_D, int T_Q>
struct CollisionParameters {
	double density = 0.0;
	double velocity[T_D] = { 0.0 };
	double scaling = 0.0;
	double tau = 0.0;
	const Stencil& stencil;
	// Unit vector for the stencil direction that consists of [1,0,-1] instead of [scaling,0,-scaling]
	double e[T_Q][T_D] = {{0.0}};
	double weight[T_Q] = {0.0};
	//unscaled speed of sound
	double cs2 = 0.0;
	double dt = 0.0;

	// Individual parameters that are needed for specific collision models only

	CollisionParameters(double scaling, double viscosity, const Stencil& st,
			double scaled_cs2, double dt) :
			scaling(scaling), stencil(st), dt(dt) {

		assert(st.getD() == T_D);
			assert(st.getQ() == T_Q);

		for (int i = 0; i < T_D; ++i) {
			for (int j = 0; j < T_Q; ++j) {
				e[j][i] = st.getDirections().at(j)(i) / scaling;
			}
		}

		for (int j = 0; j < T_Q; ++j) {
						weight[j] = st.getWeight(j);
					}

		//unscale the speed of sound
		cs2 = scaled_cs2 / (scaling * scaling);

		// The relaxation time has to be calculated with the scaled speed of sound
		tau = calculateTauFromNu(viscosity, scaled_cs2, dt);


		assert((cs2 - 1./3.) < 1e-10);

	}
};

template<int T_D, int T_Q>
inline void calculateVelocity(const double fLocal[], double velocity[],
		double scaling, double density, CollisionParameters<T_D,T_Q>& params) {
	for (int j = 0; j < T_D; j++) {
		velocity[j] = 0.0;
		for (int i = 0; i < T_Q; i++) {
			velocity[j] += params.e[i][j] * fLocal[i];
		}
		velocity[j] = velocity[j] * 1.0 / density;
	}
}

template<>
inline void calculateVelocity<2, 9>(const double fLocal[], double velocity[],
		double scaling, double density, CollisionParameters<2,9>& params) {
	velocity[0] = 1.0 / density
			* (fLocal[1] + fLocal[5] + fLocal[8] - fLocal[3] - fLocal[6]
					- fLocal[7]);
	velocity[1] = 1.0 / density
			* (fLocal[2] + fLocal[5] + fLocal[6] - fLocal[4] - fLocal[7]
					- fLocal[8]);
}
template<>
inline void calculateVelocity<3, 19>(const double fLocal[], double velocity[],
		double scaling, double density, CollisionParameters<3, 19>& params) {

	velocity[0] = 1.0 / density
			* (fLocal[1] - fLocal[3] + fLocal[7] - fLocal[8] - fLocal[9]
					+ fLocal[10] + fLocal[11] + fLocal[12] - fLocal[13]
					- fLocal[14]);
	velocity[1] = 1.0 / density
			* (-fLocal[5] + fLocal[6] - fLocal[11] + fLocal[12] + fLocal[13]
					- fLocal[14] - fLocal[15] + fLocal[16] + fLocal[17]
					- fLocal[18]);
	velocity[2] = 1.0 / density
			* (fLocal[2] - fLocal[4] + fLocal[7] + fLocal[8] - fLocal[9]
					- fLocal[10] + fLocal[15] + fLocal[16] - fLocal[17]
					- fLocal[18]);
}

} /* namespace natrium */

#endif /* AUXILIARYCOLLISIONFUNCTIONS_H_ */
