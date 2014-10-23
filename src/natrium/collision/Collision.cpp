/**
 * @file CollisionModel.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Collision.h"

#include "../solver/DistributionFunctions.h"

#include "../boltzmannmodels/D2Q9IncompressibleModel.h"
#include "BGKTransformed.h"

#include <cmath>

namespace natrium {

// constructor
template<class BoltzmannType, class CollisionType>
Collision<BoltzmannType, CollisionType>::Collision(
		BoltzmannType boltzmannType,
		CollisionType collisionType) :
		m_boltzmannType(boltzmannType), m_collisionType(collisionType), m_d(
				boltzmannType.getD()), m_q(boltzmannType.getQ()) {

} // constructor
template Collision<D2Q9IncompressibleModel, BGKTransformed>::Collision(
		D2Q9IncompressibleModel boltzmannType,
		BGKTransformed collisionModel);

template<class BoltzmannType, class CollisionType>
Collision<BoltzmannType, CollisionType>::~Collision() {
}
template Collision<D2Q9IncompressibleModel, BGKTransformed>::~Collision();

template<class BoltzmannType, class CollisionType>
void Collision<BoltzmannType, CollisionType>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = m_boltzmannType.getQ();

	assert(f.size() == Q);
	assert(velocities.size() == m_d);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

//for all degrees of freedom
	for (size_t i = 0; i < n_dofs; i++) {

		// calculate density
		densities(i) = 0;
		for (size_t j = 0; j < Q; j++) {
			densities(i) += f.at(j)(i);
		}
		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate velocity
			// for all velocity components
			for (size_t j = 0; j < m_d; j++) {
				velocities.at(j)(i) = 0;

				for (size_t k = 0; k < Q; k++) {
					velocities.at(j)(i) += f.at(k)(i)
							* m_boltzmannType.getDirection(k)(j);
				}
				velocities.at(j)(i) /= densities(i);
			}
		}

		// calculate equilibrium distribution
		// TODO Optimize by passing different arguments
		vector<double> feq(Q, 0.0);
		numeric_vector u(m_d);
		for (size_t j = 0; j < m_d; j++) {
			u(j) = velocities.at(j)(i);
		}
		m_boltzmannType.getEquilibriumDistributions(feq, u, densities(i));

		// BGK collision
		m_collisionType.collideSingleDoF(i, feq, f);
	}
}
template
void Collision<D2Q9IncompressibleModel, BGKTransformed>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const;

template<class BoltzmannType, class CollisionType>
void Collision<BoltzmannType, CollisionType>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities, const vector<bool>& isBoundary,
		bool inInitializationProcedure) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = m_boltzmannType.getQ();

	assert(f.size() == Q);
	assert(velocities.size() == m_d);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

//for all degrees of freedom
	for (size_t i = 0; i < n_dofs; i++) {

		// calculate density
		densities(i) = 0;
		for (size_t j = 0; j < Q; j++) {
			densities(i) += f.at(j)(i);
		}
		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate velocity
			// for all velocity components
			for (size_t j = 0; j < m_d; j++) {
				velocities.at(j)(i) = 0;

				for (size_t k = 0; k < Q; k++) {
					velocities.at(j)(i) += f.at(k)(i)
							* m_boltzmannType.getDirection(k)(j);
				}
				velocities.at(j)(i) /= densities(i);
			}
		}

		if (isBoundary.at(i)) {
			continue;
		}

		// calculate equilibrium distribution
		// TODO Optimize by passing different arguments
		vector<double> feq(Q, 0.0);
		numeric_vector u(m_d);
		for (size_t j = 0; j < m_d; j++) {
			u(j) = velocities.at(j)(i);
		}
		m_boltzmannType.getEquilibriumDistributions(feq, u, densities(i));

		// BGK collision
		m_collisionType.collideSingleDoF(i, feq, f);
	}
} /*collideAll */
template void Collision<D2Q9IncompressibleModel, BGKTransformed>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities, const vector<bool>& isBoundary,
		bool inInitializationProcedure) const;

} /* namespace natrium */
