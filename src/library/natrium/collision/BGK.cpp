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
		CollisionModel(stencil), m_relaxationParameter(relaxationParameter), m_prefactor(
				-1. / (relaxationParameter + 0.5)), m_dt(dt) {

} // constructor

// destructor
BGK::~BGK() {
} // destructor

void BGK::collideSinglePoint(vector<double>& distributions) const {

	// assert
	size_t Q = getStencil()->getQ();
	assert(distributions.size() == Q);

	// calculate macroscopic entities
	double rho = this->calculateDensity(distributions);
	numeric_vector u(getStencil()->getD());
	this->calculateVelocity(distributions, rho, u);

	// calculate equilibrium distribution (feq)
	vector<double> feq(Q);
	this->getEquilibriumDistributions(feq, u, rho);

	// update distribution
	for (size_t i = 0; i < Q; i++) {
		distributions.at(i) = distributions.at(i)
				+ getPrefactor() * (distributions.at(i) - feq.at(i));
	}

} //collideSinglePoint

void BGK::collideAll(DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	LOG(WARNING)
			<< "Warning: You are using a collision model that has no efficient implementation."
					"Therefore, collisions can well be as costly as streaming."
					"Providing an efficient implementation could give a speed-up of nearly 2."
			<< endl;

	size_t Q = getStencil()->getQ();
	size_t D = getStencil()->getD();

	assert(f.size() == Q);
	assert(velocities.size() == D);

#ifdef DEBUG
	size_t n_dofs = f.at(0).size();
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < D; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;

		// calculate density
		densities(i) = 0;
		for (size_t j = 0; j < Q; j++) {
			densities(i) = densities(i) + f.at(j)(i);
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
					velocities.at(j)(i) = velocities.at(j)(i)
							+ f.at(k)(i) * getStencil()->getDirection(k)(j);
				}
				velocities.at(j)(i) = velocities.at(j)(i) / densities(i);
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

