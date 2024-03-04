/*
 * AuxiliaryCollisionFunctions.h
 *
 *  Created on: 16.01.2017
 *      Author: natrium
 */
#ifndef AUXILIARYCOLLISIONFUNCTIONS_H_
#define AUXILIARYCOLLISIONFUNCTIONS_H_

#include <array>
#include <vector>

#include "../collision/ExternalForceFunctions.h"
#include "../utilities/NATriuMException.h"
#include "../utilities/ConfigNames.h"
#include "../solver/DistributionFunctions.h"
#include "../solver/SolverConfiguration.h"


namespace natrium {

// forward declaration
class CollisionException;

/**
 * @short Exception class for Unstable Collision
 */
class DensityZeroException: public NATriuMException {
private:
	std::string message;
public:
	explicit DensityZeroException(const char *msg) :
        NATriuMException(msg), message(msg) {
	}
    explicit DensityZeroException(const string& msg) :
	    NATriuMException(msg), message(msg) {
	}
	~DensityZeroException() noexcept {
	}
	const char *what() const noexcept {
		return this->message.c_str();
	}
};

template<int T_Q>
inline double calculateDensity(const std::array<double, T_Q>& fLocal) {
	double density = 0.0;
	for (size_t p = 0; p < T_Q; ++p) {
        density += fLocal[p];
    }
	if (density < 1e-10) {
		throw CollisionException("Densities too small (< 1e-10) for collisions. Decrease time step size.");
	}
	return density;
}

inline double calculateTauFromNu(double viscosity, double cs2,
		double timeStepSize) {
	double tau;
	tau = (viscosity) / (timeStepSize * cs2) + 0.5;
	return tau;
}

/**
 * @short for BGK steady state simulations
 */
inline double calculateTauFromNuAndgamma(double viscosity, double cs2,
		double timeStepSize, double gamma) {
	double tau;
	tau = (viscosity) / (gamma * timeStepSize * cs2) + 0.5;
	return tau;
}

// Apply sutherland's law of temperature dependent viscosity
inline double calculateTauFromNuAndT(double viscosity, double cs2,
        double timeStepSize, double T) {
    double tau;
    tau = (viscosity)*sqrt(T) / (timeStepSize * cs2) + 0.5;
    return tau;
}

template<int T_Q>
inline void copyGlobalToLocalF(std::array<double, T_Q>& fLocal,
		const DistributionFunctions& f, const size_t i) {
	for (size_t p = 0; p < T_Q; ++p) {
		fLocal[p] = f.at(p)(i);
	}
}

template<int T_Q>
inline void copyLocalToGlobalF(std::array<double, T_Q>& fLocal,
		DistributionFunctions& f, const size_t i) {
	for (size_t p = 0; p < T_Q; ++p) {
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
        // the local g is stored in this array
	std::array<double, T_Q> gLocal = { };
	// the local g is stored in this array
	std::array<double, T_Q> geq = { };

    std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_Q> H3;
    std::array<std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_D>,T_Q> H4;

	double density = 0.0;
	double temperature = 1.0;
    double maskShockSensor = 0.0;

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
	// preconditioning parameter for steady-state simulations (Guo 2004)
	double gamma_steadystate = 1.0;


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
		assert(scaling == st.getScaling());

		for (size_t i = 0; i < T_D; ++i) {
			for (size_t j = 0; j < T_Q; ++j) {
				e[j][i] = st.getDirections().at(j)(i) / scaling;
			}
		}

		for (size_t j = 0; j < T_Q; ++j) {
			weight[j] = st.getWeight(j);
		}

		//unscale the speed of sound
		cs2 = scaled_cs2 / (scaling * scaling);

		// The relaxation time has to be calculated with the scaled speed of sound
		// As a dimensionless parameter, it is independent of scaling
		// (Note that the scaling affects the time step size)
		tau = calculateTauFromNu(viscosity, scaled_cs2, dt);

		if ((pd.hasExternalForce()) and (forcetype != NO_FORCING)) {
			for (size_t i = 0; i < T_D; ++i) {
				forces[i] = pd.getExternalForce()->getForce()[i];
			}
		}
		//assert((cs2 - 1. / 3.) < 1e-10);

		// for simulations with steady state equilibrium (Guo 2004)
		gamma_steadystate = 1;
		if (cfg.getEquilibriumScheme() == STEADYSTATE_EQUILIBRIUM) {
			gamma_steadystate = cfg.getBGKSteadyStateGamma();
			tau = calculateTauFromNuAndgamma(viscosity, scaled_cs2, dt,
					gamma_steadystate);
			if ((pd.hasExternalForce()) and (forcetype != NO_FORCING)) {
				for (size_t i = 0; i < T_D; ++i) {
					forces[i] = pd.getExternalForce()->getForce()[i]
							/ gamma_steadystate;
				}
			}
		}

	}

};

    template<size_t T_D, size_t T_Q>
    inline std::array<std::array<double, T_D>, T_Q> getParticleVelocitiesWithoutScaling(const Stencil &st) {
        std::array<std::array<double, T_D>, T_Q> e;
        for (size_t i = 0; i < T_D; ++i) {
            for (size_t j = 0; j < T_Q; ++j) {
                e[j][i] = st.getDirections().at(j)(i) / st.getScaling();
            }
        }
        return e;
    }



template <size_t T_D>
constexpr std::array<std::array<size_t,T_D>, T_D> unity_matrix()
{
    std::array<std::array<size_t,T_D>, T_D> eye = {};
    for (size_t a = 0; a<T_D; a++)  {
        for (size_t b = 0; b < T_D; b++) {
            if (a==b){
                eye[a][b] = 1;
            }
        }
    }
    return eye;
}

template<int T_D, int T_Q>
inline void calculateVelocity(const std::array<double, T_Q>& fLocal,
		std::array<double, T_D>& velocity, double density,
		GeneralCollisionData<T_D, T_Q>& params) {
    for (size_t j = 0; j < T_D; j++) {
		velocity[j] = 0.0;
		for (size_t i = 0; i < T_Q; i++) {
			velocity[j] += params.e[i][j] * fLocal[i];
		}
		velocity[j] = velocity[j] * 1.0 / density;
	}
}

template<int T_D, int T_Q>
inline void calculateVelocity(const std::array<double, T_Q>& fLocal,
                              std::array<double, T_D>& velocity, double density,
                              std::array<std::array<double, T_D>, T_Q> e) {
    for (size_t j = 0; j < T_D; j++) {
        velocity[j] = 0.0;
        for (size_t i = 0; i < T_Q; i++) {
            velocity[j] += e[i][j] * fLocal[i];
        }
        velocity[j] = velocity[j] * 1.0 / density;
    }
}


template<>
inline void calculateVelocity<2, 9>(const std::array<double, 9>& fLocal,
		std::array<double, 2>& velocity, double density,
		GeneralCollisionData<2, 9>& params) {
    (void)params;
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
    (void)params;
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

    template<int T_D, int T_Q>
    inline double calculateTemperature(const std::array<double, T_Q> &fLocal, const std::array<double, T_Q> &gLocal,
                                       std::array<double, T_D> &velocity, double density, double temperature,
                                       GeneralCollisionData <T_D, T_Q> &params, double d) {
        (void) d;
        //T0[i,j]+=((c[k,0]-u[0,i,j])**2+(c[k,1]-u[1,i,j])**2)*fin[k,i,j]*0.5/rho[i,j]
        temperature = 0.0;
        for (size_t i = 0; i < T_Q; i++) {
            double sum = 0.0;
            for (size_t a = 0; a < T_D; a++) {
                sum += (params.e[i][a] - velocity[a]) * (params.e[i][a] - velocity[a]);
            }
            temperature += sum * fLocal[i] / params.cs2 +
                           gLocal[i];
        }
        const double gamma = params.configuration.getHeatCapacityRatioGamma();
        const double C_v = 1./(gamma-1.0);
        temperature = temperature * 0.5 / (density*C_v);
return temperature;
}

    template<int T_D, int T_Q>
    inline double calculateTemperature(const std::array<double, T_Q> &fLocal, const std::array<double, T_Q> &gLocal,
                                       const std::array<double, T_D> &velocity, const double density,
                                       const std::array<std::array<double, T_D>, T_Q> e, const double cs2, const double gamma) {

        double temperature = 0.0;
        for (size_t i = 0; i < T_Q; i++) {
            double sum = 0.0;
            for (size_t a = 0; a < T_D; a++) {
                sum += (e[i][a] - velocity[a]) * (e[i][a] - velocity[a]);
            }
            temperature += sum * fLocal[i] / cs2 +
                           gLocal[i];
        }
        const double C_v = 1./(gamma-1.0);
        temperature = temperature * 0.5 / (density*C_v);
        return temperature;
    }




template<size_t T_D, size_t T_Q>
inline void applyMacroscopicForces(std::array<double*, T_D>& velocities,
		int ii, GeneralCollisionData<T_D, T_Q>& genData) {
	assert(velocities.size() == T_D);
	if (genData.forcetype == NO_FORCING) {
		throw NATriuMException(
				"Problem requires forcing scheme, but forcing was switched off."
						"Please set forcing to SHIFTING_VELOCITY in the Solver Configuration.");
	}

	if (genData.forcetype == SHIFTING_VELOCITY) {
        // TODO: incorporate into calculate velocities  (accessing the global velocity vector twice is ugly)
        // 		 upon refactoring this, remember to incorporate the test for problem.hasExternalForce() (cf. collideAll)
        for (size_t j = 0; j < T_D; j++) {
            velocities[j][ii] = velocities[j][ii]
                                + 0.5 * genData.dt * genData.forces[j] / genData.density;
            // I wanted this to be independent of the order of execution with applyForces  (therefor 2x acces to velocities; += is risky for TrilinosVector)
        }
    }
    else if (genData.forcetype == EXACT_DIFFERENCE) {
            // PASS


	} else {
		throw NotImplementedException(
				"Force Type not implemented. Use Shifting Velocity instead.");
	}
}

template<size_t T_D, size_t T_Q>
inline void applyForces(GeneralCollisionData<T_D, T_Q>& genData) {
	if (genData.forcetype == NO_FORCING) {
		throw NATriuMException(
				"Problem requires forcing scheme, but forcing was switched off."
						"Please set forcing to SHIFTING_VELOCITY in the Solver Configuration.");
	}
	if (genData.forcetype == SHIFTING_VELOCITY) {
        for (size_t i = 0; i < T_D; i++) {
            genData.velocity[i] += genData.tau * genData.dt * genData.forces[i]
                                   / genData.density / genData.scaling;

        }
    }

    else if (genData.forcetype == EXACT_DIFFERENCE) {
        // PASS

    }

	 else {
		throw NotImplementedException(
				"Force Type not implemented. Use Shifting Velocity instead.");
	}
}

    template<int T_D, int T_Q, template<int, int> class T_equilibrium>
    inline void postCollisionApplyForces(std::array<double *, T_D> &velocities,
                                         int ii, GeneralCollisionData<T_D, T_Q> &genData) {
        if (genData.forcetype == NO_FORCING) {
            throw NATriuMException(
                    "Problem requires forcing scheme, but forcing was switched off."
                    "Please set forcing to SHIFTING_VELOCITY in the Solver Configuration.");
        }
        if (genData.forcetype == SHIFTING_VELOCITY) {
            // PASS
        } else if (genData.forcetype == EXACT_DIFFERENCE) {
            for (size_t j = 0; j < T_D; j++) {
                // no tau!
                genData.velocity[j] += genData.dt * genData.forces[j]
                                       / genData.density / genData.scaling;
                velocities[j][ii] = velocities[j][ii]
                                    + 0.5 * genData.dt * genData.forces[j] / genData.density;
            }
            std::array<double, T_Q> shiftedEq = {};
            T_equilibrium<T_D, T_Q> eq_shifted(genData);
            eq_shifted.calc(shiftedEq, genData);
            for (size_t i = 0; i < T_Q; i++) {
                genData.fLocal[i] += (shiftedEq[i] - genData.feq[i]);
            }


        } else {
            throw NotImplementedException(
                    "Force Type not implemented. Use Shifting Velocity instead.");
        }

    }


template<size_t T_D, size_t T_Q>
inline void calculateGeqFromFeq(const std::array<double, T_Q>& feq,std::array<double, T_Q>& geq, const GeneralCollisionData<T_D,T_Q>& genData)
{
    const double gamma = genData.configuration.getHeatCapacityRatioGamma();
    const double C_v = 1. / (gamma - 1.0);
    for (size_t i = 0; i < T_Q; i++) {
        geq[i]=feq[i]*(genData.temperature)*(2.0*C_v-T_D);
    }
}

template<size_t T_D, size_t T_Q>
inline void calculateGeqFromFeq(const std::array<double, T_Q>& feq,std::array<double, T_Q>& geq, const double temperature, const double gamma)
{
    const double C_v = 1. / (gamma - 1.0);
    for (size_t i = 0; i < T_Q; i++) {

        geq[i]=feq[i]*(temperature)*(2.0*C_v-T_D);

    }

}

    template<size_t T_D, size_t T_Q>
    inline void calculateCenteredHeatFluxTensor(const std::array<double, T_Q> &f,
                                                std::array<std::array<std::array<double, T_D>, T_D>, T_D> &heatFluxTensor,
                                                const GeneralCollisionData<T_D, T_Q> &p) {
        for (size_t i = 0; i < T_Q; i++) {
            for (size_t a = 0; a < T_D; a++) {
                for (size_t b = 0; b < T_D; b++) {
                    for (size_t c = 0; c < T_D; c++) {

                        heatFluxTensor[a][b][c] += ((p.e[i][a] - p.velocity[a]) * (p.e[i][b] - p.velocity[b]) *
                                                    (p.e[i][c] - p.velocity[c])) * f[i];
                    }
                }
            }
        }
    }

    template<size_t T_D, size_t T_Q>
    inline void calculateCenteredMomentumFlux(const std::array<double, T_Q> &gNeq,
                                                std::array<double, T_D> &CentralFlux,
                                                const GeneralCollisionData<T_D, T_Q> &p) {

        for (size_t i = 0; i < T_Q; i++) {
            for (size_t a = 0; a < T_D; a++) {
                CentralFlux[a] += (p.e[i][a] - p.velocity[a]) * gNeq[i] ;
            }
        }
    }

    template<size_t T_D, size_t T_Q>
    inline double calculateKnudsenNumberEstimate(const std::array<double, T_Q> &f, const std::array<double, T_Q> &feq, const std::array<double, T_Q> &weight) {
    double estimate = 0.0;
        for (size_t i = 0; i < T_Q; i++) {
            estimate += abs(f[i] - feq[i]) / weight[i];
    }
    return estimate / T_Q;
    }


    template<size_t T_D, size_t T_Q>
    inline void
    calculateFStar(std::array<double, T_Q> &fStar, const std::array<std::array<std::array<double, T_D>, T_D>, T_D> &QNeq,
                   const GeneralCollisionData<T_D, T_Q> &p) {
        std::array<std::array<size_t,T_D>, T_D> eye = unity_matrix<T_D>();
        const double cs6 = 6.0 * p.cs2 * p.cs2 * p.cs2;
            for (size_t a = 0; a < T_D; a++) {
                for (size_t b = 0; b < T_D; b++) {
                    for (size_t c = 0; c < T_D; c++) {
                        for (size_t i = 0; i < T_Q; i++) {
                        fStar[i] += p.weight[i] * (QNeq[a][b][c] * (p.e[i][a] * p.e[i][b] * p.e[i][c] -
                                                                    3 * p.cs2 * p.e[i][c] * eye[a][b])) /
                                    cs6;
                    }
                }
            }
        }
    }

    template<size_t T_D, size_t T_Q>
    inline void
    calculateGStar(std::array<double, T_Q> &gStar, const std::array<double, T_D> &centeredFluxTensorG,
                   const GeneralCollisionData<T_D, T_Q> &p) {
//        std::array<std::array<size_t,T_D>, T_D> eye = unity_matrix<T_D>();
//        const double cs6 = 6.0 * p.cs2 * p.cs2 * p.cs2;
        for (size_t a = 0; a < T_D; a++) {
            for (size_t i = 0; i < T_Q; i++) {
                gStar[i] += p.weight[i] * (centeredFluxTensorG[a] * p.e[i][a]) /
                                    p.temperature;
                    }
                }
            }

    template<size_t T_D, size_t T_Q>
    inline std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_Q> calculateH3(const double cs2, std::array<std::array<double,T_D>,T_Q> e) {
        std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_Q> H3;
        const std::array<std::array<size_t,T_D>, T_D> eye = unity_matrix<T_D>();

        for (size_t i = 0; i < T_Q; i++) {
            for (size_t a = 0; a < T_D; a++) {
                for (size_t b = 0; b < T_D; b++) {
                    for (size_t c = 0; c < T_D; c++) {
                        H3[i][a][b][c] = e[i][a] * e[i][b] * e[i][c]
                                            - cs2 * (e[i][a]*eye[b][c] + e[i][b]*eye[a][c] + e[i][c]*eye[a][b]);
                    }
                }
            }
        }
        return H3;

    }

    template<size_t T_D, size_t T_Q>
    inline std::array<std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_D>,T_Q> calculateH4(const double cs2, std::array<std::array<double,T_D>,T_Q> e) {
        std::array<std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_D>,T_Q> H4;
        const std::array<std::array<size_t,T_D>, T_D> eye = unity_matrix<T_D>();

        for (size_t i = 0; i < T_Q; i++) {
            for (size_t a = 0; a < T_D; a++) {
                for (size_t b = 0; b < T_D; b++) {
                    for (size_t c = 0; c < T_D; c++) {
                        for (size_t d = 0; d < T_D; d++) {

                            const double power4 = e[i][a] * e[i][b] * e[i][c] * e[i][d];
                            const double power2 = e[i][a] * e[i][b] * eye[c][d]
                                            + e[i][a] * e[i][c] * eye[b][d]
                                            + e[i][a] * e[i][d] * eye[b][c]
                                            + e[i][b] * e[i][c] * eye[a][d]
                                            + e[i][b] * e[i][d] * eye[a][c]
                                            + e[i][c] * e[i][d] * eye[a][b];
                            const double power0 = eye[a][b] * eye[c][d] + eye[a][c] * eye[b][d] +
                                            eye[a][d] * eye[b][c];

                            H4[i][a][b][c][d] = (power4 - cs2 * power2 + cs2 * cs2 * power0);
                        }
                    }
                }
            }
        }
        return H4;

    }

} /* namespace natrium */

#endif /* AUXILIARYCOLLISIONFUNCTIONS_H_ */
