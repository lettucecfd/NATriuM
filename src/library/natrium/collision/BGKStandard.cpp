/**
 * @file BGKStandard.cpp
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGKStandard.h"

namespace natrium {

/// constructor
BGKStandard::BGKStandard(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		BGK(relaxationParameter, dt, stencil) {
} /// constructor

/// destructor
BGKStandard::~BGKStandard() {
} /// destructor

double BGKStandard::getEquilibriumDistribution(size_t i,
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

void BGKStandard::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else if (Stencil_D3Q19 == getStencil()->getStencilType()) {
		collideAllD3Q19(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else {
		// Inefficient collision
		BGK::collideAll(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	}
}

void BGKStandard::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	// Efficient collision for D2Q9
	size_t Q = 9;
	size_t D = 2;
	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare();
	double prefactor = scaling / cs2;
	double relax_factor = getPrefactor();

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

	// allocation
	double feq[9] = { };
	double scalar_product;
	double uSquareTerm;
	double mixedTerm;
	double weighting;
	double rho_i;
	double u_0_i;
	double u_1_i;
	double f_i[9] = { };

	distributed_vector& f0 = f.at(0);
	distributed_vector& f1 = f.at(1);
	distributed_vector& f2 = f.at(2);
	distributed_vector& f3 = f.at(3);
	distributed_vector& f4 = f.at(4);
	distributed_vector& f5 = f.at(5);
	distributed_vector& f6 = f.at(6);
	distributed_vector& f7 = f.at(7);
	distributed_vector& f8 = f.at(8);

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// copy fi
		f_i[0] = f0(i);
		f_i[1] = f1(i);
		f_i[2] = f2(i);
		f_i[3] = f3(i);
		f_i[4] = f4(i);
		f_i[5] = f5(i);
		f_i[6] = f6(i);
		f_i[7] = f7(i);
		f_i[8] = f8(i);

		// calculate density
		rho_i = f_i[0] + f_i[1] + f_i[2] + f_i[3]
				+ f_i[4] + f_i[5] + f_i[6] + f_i[7]
				+ f_i[8];
		densities(i) = rho_i;
		if (rho_i < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate velocity
			// for all velocity components
			u_0_i = scaling / rho_i
					* (f_i[1] + f_i[5] + f_i[8] - f_i[2]
							- f_i[6] - f_i[7]);
			u_1_i = scaling / rho_i
					* (f_i[2] + f_i[5] + f_i[6] - f_i[4]
							- f_i[7] - f_i[8]);
			velocities.at(0)(i) = u_0_i;
			velocities.at(1)(i) = u_1_i;
		}

		// calculate equilibrium distribution
		scalar_product = u_0_i * u_0_i
				+ u_1_i * u_1_i;
		uSquareTerm = -scalar_product / (2 * cs2);
		// direction 0
		weighting = 4. / 9. * rho_i;
		feq[0] = weighting * (1 + uSquareTerm);
		// directions 1-4
		weighting = 1. / 9. * rho_i;
		mixedTerm = prefactor * (u_0_i);
		feq[1] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[3] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (u_1_i);
		feq[2] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[4] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		// directions 5-8
		weighting = 1. / 36. * rho_i;
		mixedTerm = prefactor * (u_0_i + u_1_i);
		feq[5] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[7] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (-u_0_i + u_1_i);
		feq[6] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[8] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);

		// BGK collision
		f0(i) += relax_factor * (f_i[0] - feq[0]);
		f1(i) += relax_factor * (f_i[1] - feq[1]);
		f2(i) += relax_factor * (f_i[2] - feq[2]);
		f3(i) += relax_factor * (f_i[3] - feq[3]);
		f4(i) += relax_factor * (f_i[4] - feq[4]);
		f5(i) += relax_factor * (f_i[5] - feq[5]);
		f6(i) += relax_factor * (f_i[6] - feq[6]);
		f7(i) += relax_factor * (f_i[7] - feq[7]);
		f8(i) += relax_factor * (f_i[8] - feq[8]);
	} /* for all dofs */
} /* collideAllD2Q9 */

void BGKStandard::collideAllD3Q19(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

} /* collideAllD3Q19 */

void BGKStandard::collideAllD3Q15(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {
	assert (false);
}

} /* namespace natrium */
