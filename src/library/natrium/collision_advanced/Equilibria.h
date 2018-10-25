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
#include <array>


namespace natrium{
template <int T_D,int T_Q>
class BGKEquilibrium
{

public:
	//BGKEquilibrium();
	void calc(std::array<double, T_Q>& feq, const GeneralCollisionData<T_D,T_Q> & params);
};



template<>
 inline void BGKEquilibrium<2,9>::calc(std::array<double, 9>& feq, const GeneralCollisionData<2,9> & params) {
	double scalar_product;
	double weighting;
	double uSquareTerm;
	double mixedTerm;
	double prefactor = 1. / params.cs2;
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

template<>
 inline void BGKEquilibrium<2,25>::calc(std::array<double, 25>& feq, const GeneralCollisionData<2,25> & params)
	{
double eye[2][2]={{1,0},{0,1}};
		double uu_term = 0.0;
		for (size_t j = 0; j < 2; j++) {
			uu_term += -(params.velocity[j] * params.velocity[j])
					/ (2.0 * params.cs2);
		}

		for (size_t i = 0; i < 25; i++) {
			double ue_term = 0.0;
			for (size_t j = 0; j < 2; j++) {
				ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
			}
			feq[i] = params.weight[i] * params.density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
			for (int alp = 0; alp < 2; alp++){
				for (int bet = 0; bet < 2; bet++){
					for (int gam = 0; gam < 2; gam++){
						feq[i] += params.weight[i] * params.density * params.velocity[alp]*params.velocity[bet]*params.velocity[gam]*params.e[i][gam]/(6.*params.cs2*params.cs2*params.cs2)*(params.e[i][alp]*params.e[i][bet]-3.*params.cs2*eye[alp][bet]);
					}
				}

			}
		}
	}



template<int T_D, int T_Q>
inline void BGKEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q>& feq,
		const GeneralCollisionData<T_D, T_Q> & params) {

	double uu_term = 0.0;
	for (size_t j = 0; j < T_D; j++) {
		uu_term += -(params.velocity[j] * params.velocity[j])
				/ (2.0 * params.cs2);
	}

	for (size_t i = 0; i < T_Q; i++) {
		double ue_term = 0.0;
		for (size_t j = 0; j < T_D; j++) {
			ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
		}
		feq[i] = params.weight[i] * params.density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
	}
} /* BGKEquilibrium<T_D, T_Q>::calc */

} /* class BGKEquilibrium */




// ==========================================================================================

namespace natrium{
template <int T_D,int T_Q>
class SteadyStateEquilibrium
{

public:
	void calc(std::array<double, T_Q>& feq, const GeneralCollisionData<T_D,T_Q> & params);
};

template<int T_D, int T_Q>
inline void SteadyStateEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q>& feq,
		const GeneralCollisionData<T_D, T_Q> & params) {

	double uu_term = 0.0;
	for (size_t j = 0; j < T_D; j++) {
		uu_term += -(params.velocity[j] * params.velocity[j])
				/ (2.0 * params.cs2);
	}

	for (size_t i = 0; i < T_Q; i++) {
		double ue_term = 0.0;
		for (size_t j = 0; j < T_D; j++) {
			ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
		}
		feq[i] = params.weight[i] * params.density * (1 + ue_term * (1 + 0.5 * (ue_term) / params.gamma_steadystate)
				+ uu_term / params.gamma_steadystate);
	}
} /* SteadyStateEquilibrium<T_D, T_Q>::calc */

} /* class SteadyStateEquilibrium */

// ==========================================================================================

namespace natrium{
template <int T_D,int T_Q>
class IncompressibleEquilibrium
{

public:
	void calc(std::array<double, T_Q>& feq, const GeneralCollisionData<T_D,T_Q> & params);
};

template<int T_D, int T_Q>
inline void IncompressibleEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q>& feq,
		const GeneralCollisionData<T_D, T_Q> & params) {

	double uu_term = 0.0;
	for (size_t j = 0; j < T_D; j++) {
		uu_term += -(params.velocity[j] * params.velocity[j])
				/ (2.0 * params.cs2);
	}

	for (size_t i = 0; i < T_Q; i++) {
		double ue_term = 0.0;
		for (size_t j = 0; j < T_D; j++) {
			ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
		}
		feq[i] = params.weight[i] * (params.density + ue_term * (1 + 0.5 * (ue_term) )
				+ uu_term );
	}
} /* IncompressibleEquilibrium<T_D, T_Q>::calc */

} /* class IncompressibleEquilibrium */


// ==========================================================================================

namespace natrium{
template <int T_D,int T_Q>
class EntropicEquilibrium
{

public:
	void calc(std::array<double, T_Q>& feq, const GeneralCollisionData<T_D,T_Q> & params);
};

template<int T_D, int T_Q>
inline void EntropicEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q>& feq,
		const GeneralCollisionData<T_D, T_Q> & params) {

	throw NotImplementedException("Entropic equilibrium is not implemented, yet.");
	/*double uu_term = 0.0;
	for (size_t j = 0; j < T_D; j++) {
		uu_term += -(params.velocity[j] * params.velocity[j])
				/ (2.0 * params.cs2);
	}

	for (size_t i = 0; i < T_Q; i++) {
		double ue_term = 0.0;
		for (size_t j = 0; j < T_D; j++) {
			ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
		}
		feq[i] = params.weight[i] * (params.density + ue_term * (1 + 0.5 * (ue_term) )
				+ uu_term );
	}*/
} /* EntropicEquilibrium<T_D, T_Q>::calc */

} /* class EntropicEquilibrium */


#endif /* LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_ */
