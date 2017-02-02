/*
 * AuxiliaryCollisionFunctions.h
 *
 *  Created on: 16.01.2017
 *      Author: natrium
 */
#ifndef AUXILIARYCOLLISIONFUNCTIONS_H_
#define AUXILIARYCOLLISIONFUNCTIONS_H_

#include "../collision/ExternalForceFunctions.h"

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
inline void copyLocalToGlobalF(double fLocal[], const DistributionFunctions& f,
		const size_t i) {
	for (int p = 0; p < T_Q; ++p) {
	//	f[p][i] = fLocal[p];
	}
}
template<int T_D, int T_Q>
struct CollisionParameters {
	double density = 0.0;
	double velocity[T_D] = { 0.0 };
	double scaling = 0.0;
	double tau = 0.0;
	double cs2 = 0.0;
	double dt = 0.0;
	CollisionParameters(double scaling, double viscosity, double cs2, double dt) :
			scaling(scaling), cs2(cs2), dt(dt) {
		tau = calculateTauFromNu(viscosity, cs2, dt);
	}
	CollisionParameters(){}
};

template<int T_D, int T_Q>
inline void calculateVelocity(const double fLocal[], double velocity[],
		double scaling, double density);

template<>
inline void calculateVelocity<2, 9>(const double fLocal[], double velocity[],
		double scaling, double density) {
	velocity[0] = scaling / density
			* (fLocal[1] + fLocal[5] + fLocal[8] - fLocal[3] - fLocal[6]
					- fLocal[7]);
	velocity[1] = scaling / density
			* (fLocal[2] + fLocal[5] + fLocal[6] - fLocal[4] - fLocal[7]
					- fLocal[8]);
}

} /* namespace natrium */

#endif /* AUXILIARYCOLLISIONFUNCTIONS_H_ */
