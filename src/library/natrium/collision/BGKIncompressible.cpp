/**
 * @file BGKIncompressible.cpp
 * @short D2Q9 model description for improved incompressible flow by Guo et al.(2000)s.
 * @date 19.07.2015
 * @author Dominik Wilde, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGKIncompressible.h"

namespace natrium {

BGKIncompressible::BGKIncompressible(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		BGK(relaxationParameter, dt, stencil) {
}

BGKIncompressible::~BGKIncompressible() {
}

double BGKIncompressible::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {
	assert(i < getStencil()->getQ());
	assert(rho > 0);
	assert(u.size() == getStencil()->getD());
	assert(u(0) < 1000000000000000.);
	assert(u(1) < 1000000000000000.);

	double prefactor = getStencil()->getWeight(i) * rho;
	double uSquareTerm = -Math::scalar_product(u, u)
			/ (2 * getStencil()->getSpeedOfSoundSquare());
	if (0 == i) {
		return prefactor * (1 + uSquareTerm);
	}
	double mixedTerm = Math::scalar_product(u, getStencil()->getDirection(i))
			/ getStencil()->getSpeedOfSoundSquare();
	return prefactor * (1 + mixedTerm * (1 + 0.5 * (mixedTerm)) + uSquareTerm);

} /// getEquilibriumDistribution

void BGKIncompressible::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} /*else if (Stencil_D3Q19 == getStencil()->getStencilType()) {
		collideAllD3Q19(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else if (Stencil_D3Q15 == getStencil()->getStencilType()) {
		collideAllD3Q15(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} */else {
		throw CollisionException("BGKIncompressible only implemented for D2Q9");
		// Inefficient collision
		//BGK::collideAll(f, densities, velocities, locally_owned_dofs,
		//		inInitializationProcedure);
	}
}

void BGKIncompressible::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {


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
		double weighting;

		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;


			// calculate density
			densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);

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

			//BGK Incompressible parameters
			double sigma = m_lambda + m_gamma;

			// calculate equilibrium distribution
			scalar_product = velocities.at(0)(i) * velocities.at(0)(i)
					+ velocities.at(1)(i) * velocities.at(1)(i);
			uSquareTerm = -scalar_product / (2 * cs2);
			// direction 0
			weighting = 4. / 9.;
			feq.at(0) = -4 * sigma * densities(i)
					+ weighting * (1 + uSquareTerm);
			// directions 1-4
			weighting = 1. / 9.;
			mixedTerm = prefactor * (velocities.at(0)(i));
			feq.at(1) = m_lambda * densities(i)
					+ weighting
							* (1 + mixedTerm * (1 + 0.5 * mixedTerm)
									+ uSquareTerm);
			feq.at(3) = m_lambda * densities(i)
					+ weighting
							* (1 - mixedTerm * (1 - 0.5 * mixedTerm)
									+ uSquareTerm);
			mixedTerm = prefactor * (velocities.at(1)(i));
			feq.at(2) = m_lambda * densities(i)
					+ weighting
							* (1 + mixedTerm * (1 + 0.5 * mixedTerm)
									+ uSquareTerm);
			feq.at(4) = m_lambda * densities(i)
					+ weighting
							* (1 - mixedTerm * (1 - 0.5 * mixedTerm)
									+ uSquareTerm);
			// directions 5-8
			weighting = 1. / 36. * densities(i);
			mixedTerm = prefactor * (velocities.at(0)(i) + velocities.at(1)(i));
			feq.at(5) = m_gamma * densities(i)
					+ weighting
							* (1 + mixedTerm * (1 + 0.5 * mixedTerm)
									+ uSquareTerm);
			feq.at(7) = m_gamma * densities(i)
					+ weighting
							* (1 - mixedTerm * (1 - 0.5 * mixedTerm)
									+ uSquareTerm);
			mixedTerm = prefactor
					* (-velocities.at(0)(i) + velocities.at(1)(i));
			feq.at(6) = m_gamma * densities(i)
					+ weighting
							* (1 + mixedTerm * (1 + 0.5 * mixedTerm)
									+ uSquareTerm);
			feq.at(8) = m_gamma * densities(i)
					+ weighting
							* (1 - mixedTerm * (1 - 0.5 * mixedTerm)
									+ uSquareTerm);

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
} /* namespace natrium */
