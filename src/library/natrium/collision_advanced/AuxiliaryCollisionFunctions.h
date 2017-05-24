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

/**
 * @short Exception class for Unstable Collision
 */
class DensityZeroException: public NATriuMException {
private:
	std::string message;
public:
	DensityZeroException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	DensityZeroException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~DensityZeroException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

template<int T_Q>
inline double calculateDensity(const std::array<double, T_Q>& fLocal) {
	double density = 0;
	for (int p = 0; p < T_Q; ++p) {
		density += fLocal[p];
	}

	if (density < 1e-10) {
		throw CollisionException(
				"Densities too small (< 1e-10) for collisions. Decrease time step size.");
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
inline void copyGlobalToLocalF(std::array<double, T_Q>& fLocal,
		const DistributionFunctions& f, const size_t i) {
	for (int p = 0; p < T_Q; ++p) {
		fLocal[p] = f.at(p)(i);
	}
}

template<int T_Q>
inline void copyLocalToGlobalF(std::array<double, T_Q>& fLocal,
		DistributionFunctions& f, const size_t i) {
	for (int p = 0; p < T_Q; ++p) {
		distributed_vector& f0 = f.at(p);
		f0(i) = fLocal[p];
	}
}
//Stores the needed parameters for the collision phase
template<int T_D, int T_Q>
struct GeneralCollisionData {
	SolverConfiguration& configuration;
	const ProblemDescription<T_D> & problemDescription;

	// the local f is stored in this array
	std::array<double, T_Q> fLocal = { };
	// the local f is stored in this array
	std::array<double, T_Q> feq = { };

	double density = 0.0;
	std::array<double, T_D> velocity = { };
	//scaling of the calculation. All parameters are unscaled during the calculation. The macroscopic velocity has to be scaled at the end of the collision step
	double scaling = 0.0;
	//relaxation time tau
	double tau = 0.0;
	//unscaled speed of sound
	double cs2 = 0.0;
	// stencil
	const Stencil& stencil;
	// time step size
	double dt = 0.0;

	ForceType forcetype;
	std::array<double, T_D> forces;

	// Unit vector for the stencil direction that consists of [1,0,-1] instead of [scaling,0,-scaling]
	std::array<std::array<double, T_D>, T_Q> e = { { } };
	// Weights of the given stencil
	std::array<double, T_Q> weight = { };

	// Individual parameters that are needed for specific collision models only

	GeneralCollisionData(SolverConfiguration& cfg,
			const ProblemDescription<T_D>& pd, double scl, double viscosity,
			const Stencil& st, double scaled_cs2, double dt) :
			configuration(cfg), problemDescription(pd), scaling(scl), stencil(
					st), dt(dt), forcetype(cfg.getForcingScheme()) {

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

		if ((pd.hasExternalForce()) and (forcetype != NO_FORCING)) {
			for (int i = 0; i < T_D; ++i) {
				forces[i] = pd.getExternalForce()->getForce()[i] / scaling / scaling;
			}
		}
		assert((cs2 - 1. / 3.) < 1e-10);

	}
};

template<int T_D, int T_Q>
inline void calculateVelocity(const std::array<double, T_Q>& fLocal,
		std::array<double, T_D>& velocity, double density,
		GeneralCollisionData<T_D, T_Q>& params) {
	for (int j = 0; j < T_D; j++) {
		velocity[j] = 0.0;
		for (int i = 0; i < T_Q; i++) {
			velocity[j] += params.e[i][j] * fLocal[i];
		}
		velocity[j] = velocity[j] * 1.0 / density;
	}
}

template<>
inline void calculateVelocity<2, 9>(const std::array<double, 9>& fLocal,
		std::array<double, 2>& velocity, double density,
		GeneralCollisionData<2, 9>& params) {
	velocity[0] = 1.0 / density
			* (fLocal[1] + fLocal[5] + fLocal[8] - fLocal[3] - fLocal[6]
					- fLocal[7]);
	velocity[1] = 1.0 / density
			* (fLocal[2] + fLocal[5] + fLocal[6] - fLocal[4] - fLocal[7]
					- fLocal[8]);
}

template<>
inline void calculateVelocity<3, 19>(const std::array<double, 19>& fLocal,
		std::array<double, 3>& velocity, double density,
		GeneralCollisionData<3, 19>& params) {

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

template<size_t T_D, size_t T_Q>
inline void applyMacroscopicForces(vector<distributed_vector>& velocities,
		size_t i, GeneralCollisionData<T_D, T_Q> genData) {
	assert(velocities.size() == T_D);
	if (genData.forcetype == NO_FORCING) {
		throw NATriuMException(
				"Problem requires forcing scheme, but forcing was switched off."
						"Please set forcing to SHIFTING_VELOCITY in the Solver Configuration.");
	}
	if (genData.forcetype == SHIFTING_VELOCITY) {
		// TODO: incorporate into calculate velocities  (accessing the global velocity vector twice is ugly)
		// 		 upon refactoring this, remember to incorporate the test for problem.hasExternalForce() (cf. collideAll)
		for (int j = 0; j < T_D; j++) {
			velocities[j](i) = velocities[j](i) + (genData.scaling * 0.5 * genData.dt * genData.forces[j]
									/ genData.density);
		}
	} else {
		throw NotImplementedException(
				"Force Type not implemented. Use Shifting Velocity instead.");
	}
}

template<size_t T_D, size_t T_Q>
inline void applyForces(GeneralCollisionData<T_D, T_Q> genData) {
	if (genData.forcetype == NO_FORCING) {
		throw NATriuMException(
				"Problem requires forcing scheme, but forcing was switched off."
						"Please set forcing to SHIFTING_VELOCITY in the Solver Configuration.");
	}
	if (genData.forcetype == SHIFTING_VELOCITY) {
		for (int i = 0; i < T_D; i++) {
			//	genData.velocity[i] *=genData.scaling;
			// TODO: Something to look at: the '-=' was a '+=' before, but BGKStandard has '-='
			genData.velocity[i] -= genData.tau * genData.dt * genData.forces[i]
					/ genData.density;
			//	genData.velocity[i] /=genData.scaling;

		}
	} else {
		throw NotImplementedException(
				"Force Type not implemented. Use Shifting Velocity instead.");
	}
}

} /* namespace natrium */

#endif /* AUXILIARYCOLLISIONFUNCTIONS_H_ */
