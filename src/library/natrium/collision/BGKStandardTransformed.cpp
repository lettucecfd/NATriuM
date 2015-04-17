/**
 * @file BGKStandardTransformed.cpp
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGKStandardTransformed.h"

namespace natrium {

/// constructor
BGKStandardTransformed::BGKStandardTransformed(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		BGK(relaxationParameter, dt, stencil),
		m_rho0 (1.0) {
} /// constructor

/// destructor
BGKStandardTransformed::~BGKStandardTransformed() {
} /// destructor

double BGKStandardTransformed::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

	assert(i < getStencil()->getQ());
	assert(rho > 0);
	assert(i >= 0);
	assert(u.size() == getStencil()->getD());
	assert(u(0) < 1000000000000000.);
	assert(u(1) < 1000000000000000.);

	double wi =  getStencil()->getWeight(i) ;
	double prefactor = wi * rho;
	double uSquareTerm = -Math::scalar_product(u, u)
			/ (2 * getStencil()->getSpeedOfSoundSquare());
	if (0 == i) {
		return prefactor * (1 + uSquareTerm) - wi * m_rho0;
	}
	double mixedTerm = Math::scalar_product(u, getStencil()->getDirection(i))
			/ getStencil()->getSpeedOfSoundSquare();
	return prefactor * (1 + mixedTerm * (1 + 0.5 * (mixedTerm)) + uSquareTerm) - wi * m_rho0;

} /// getEquilibriumDistribution


void BGKStandardTransformed::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 != getStencil()->getStencilType()) {
		// Inefficient collision for other than D2Q9
		BGK::collideAll(f, densities, velocities, inInitializationProcedure);
	} else {
		// Efficient collision for D2Q9
		size_t n_dofs = f.at(0).size();
		size_t Q = 9;
		size_t D = 2;
		double scaling = getStencil()->getScaling();
		double cs2 = getStencil()->getSpeedOfSoundSquare();
		double prefactor = scaling / cs2;
		double relax_factor = getPrefactor();

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

		// allocation
		vector<double> feq(Q, 0.0);
		double scalar_product;
		double uSquareTerm;
		double mixedTerm;
		double wi;
		double wi_rho0;
		double weighting;

		// for all dofs
		for (size_t i = 0; i < n_dofs; i++) {

			// calculate density
			densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i) + m_rho0;
			if (densities(i) < 1e-10) {
				throw CollisionException(
						"Densities too small (< 1e-10) for collisions. Decrease time step size.");
			}

			if (not inInitializationProcedure) {
				// calculate velocity
				// for all velocity components
				velocities.at(0)(i) = scaling / densities(i)
						* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
								- f.at(6)(i) - f.at(7)(i));
				velocities.at(1)(i) = scaling / densities(i)
						* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
								- f.at(7)(i) - f.at(8)(i));
			}

			// calculate equilibrium distribution
			scalar_product = velocities.at(0)(i) * velocities.at(0)(i)
					+ velocities.at(1)(i) * velocities.at(1)(i);
			uSquareTerm = -scalar_product / (2 * cs2);
			// direction 0
			wi = 4./9.;
			weighting = wi * densities(i);
			feq.at(0) = weighting * (1 + uSquareTerm) - wi * m_rho0;
			// directions 1-4
			wi = 1. / 9.;
			weighting = wi * densities(i);
			wi_rho0 = wi * m_rho0;
			mixedTerm = prefactor * (velocities.at(0)(i));
			feq.at(1) = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			feq.at(3) = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			mixedTerm = prefactor * (velocities.at(1)(i));
			feq.at(2) = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			feq.at(4) = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			// directions 5-8
			wi = 1. / 36.;
			weighting = wi * densities(i);
			wi_rho0 = wi * m_rho0;
			mixedTerm = prefactor * (velocities.at(0)(i) + velocities.at(1)(i));
			feq.at(5) = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			feq.at(7) = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			mixedTerm = prefactor
					* (-velocities.at(0)(i) + velocities.at(1)(i));
			feq.at(6) = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;
			feq.at(8) = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm) - wi_rho0;

			// BGK collision
			f.at(0)(i) += relax_factor * (f.at(0)(i) - feq.at(0));
			f.at(1)(i) += relax_factor * (f.at(1)(i) - feq.at(1));
			f.at(2)(i) += relax_factor * (f.at(2)(i) - feq.at(2));
			f.at(3)(i) += relax_factor * (f.at(3)(i) - feq.at(3));
			f.at(4)(i) += relax_factor * (f.at(4)(i) - feq.at(4));
			f.at(5)(i) += relax_factor * (f.at(5)(i) - feq.at(5));
			f.at(6)(i) += relax_factor * (f.at(6)(i) - feq.at(6));
			f.at(7)(i) += relax_factor * (f.at(7)(i) - feq.at(7));
			f.at(8)(i) += relax_factor * (f.at(8)(i) - feq.at(8));
		}
	}
}

} /* namespace natrium */
