/**
 * @file BGK.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGK.h"

namespace natrium {

/// constructor
BGK::BGK(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		CollisionModel(stencil), m_relaxationParameter(
				relaxationParameter), m_prefactor(
				-1. / (relaxationParameter + 0.5)), m_dt(dt) {

} // constructor

// destructor
BGK::~BGK() {
} // destructor


void BGK::collideSinglePoint(vector<double>& distributions) const {

	// assert
	size_t Q = getStencil()->getQ();
	assert(distributions.size() == Q);

	// create DistributionFunctions
	DistributionFunctions f;
	f.reinit(Q, 1);
	for (size_t i = 0; i < Q; i++) {
		f.at(i)(0) = distributions.at(i);
	}

	// calculate macroscopic entities
	double rho = this->calculateDensity(distributions);
	numeric_vector u(getStencil()->getD());
	this->calculateVelocity(distributions, rho, u);

	// calculate equilibrium distribution (feq)
	vector<double> feq(Q);
	this->getEquilibriumDistributions(feq, u, rho);

	// update distribution
	collideSingleDoF(0, feq, f);
	for (size_t i = 0; i < Q; i++) {
		distributions.at(i) += getPrefactor()
				* (distributions.at(i) - feq.at(i));
	}

} //collideSinglePoint

 void BGK::collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = getStencil()->getQ();
	size_t D = getStencil()->getD();

	assert(f.size() == Q);
	assert(velocities.size() == D);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < D; i++) {
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
			for (size_t j = 0; j < D; j++) {
				velocities.at(j)(i) = 0;

				for (size_t k = 0; k < Q; k++) {
					velocities.at(j)(i) += f.at(k)(i)
							* getStencil()->getDirection(k)(j);
				}
				velocities.at(j)(i) /= densities(i);
			}
		}

		// calculate equilibrium distribution
		// TODO Optimize by passing different arguments
		vector<double> feq(Q, 0.0);
		numeric_vector u(D);
		for (size_t j = 0; j < D; j++) {
			u(j) = velocities.at(j)(i);
		}
		this->getEquilibriumDistributions(feq, u, densities(i));

		// BGK collision
		collideSingleDoF(i, feq, f);
	}
}



}
/* namespace natrium */
