/*
 * EntropicStabilized.cpp
 *
 *  Created on: 11.05.2017
 *      Author: akraem3m
 */

#include "EntropicStabilized.h"
#include <sstream>
#include <functional> // for reference wrapper
#include "../dataprocessors/PseudoEntropicStabilizer.h"

namespace natrium {

EntropicStabilized::EntropicStabilized(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil) {

}

EntropicStabilized::~EntropicStabilized() {
}

// ======================================================================================================

// premature object to allow visualization of the parameter
boost::shared_ptr<distributed_vector> EntropicStabilized::OMEGA2 =
		boost::shared_ptr<distributed_vector>();

void EntropicStabilized::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	const size_t Q = MRT::getQ();

	// initialize omega2 vector
	if (not OMEGA2) {
		OMEGA2 = boost::make_shared<distributed_vector>();
	}
	if (OMEGA2->locally_owned_elements().size() == 0) {
		OMEGA2->reinit(densities);
	}

	if (9 == Q) {
		collide<2, 9>(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else if (19 == Q) {
		collide<3, 19>(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else if (27 == Q) {
		collide<3, 27>(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else {
		std::stringstream s;
		s << "Collision model not defined for Q=" << Q << endl;
		throw CollisionException(s.str());
	}

}

template<size_t D, size_t Q>
void EntropicStabilized::collide(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	// assertion
	assert (Q == f.size());
	assert (D == velocities.size());

	// Data for collision
	EStCollisionData<D, Q> _(*getStencil(), getRelaxationParameter() + 0.5);

	// store references to distribution functions to allow quick access
	std::vector<std::reference_wrapper<distributed_vector> > F;
	for (size_t i = 0; i < Q; i++) {
		F.push_back(f.at(i));
	}

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// copy f_i
		for (size_t j = 0; j < Q; j++) {
			_.f_i[j] = F[j](i);
		}

		// calculate density, momentum, velocity
		_.rho_i = density<Q>(_.f_i);
		momentum<D, Q>(_.f_i, _.j_i, _.e);
		for (size_t j = 0; j < D; j++) {
			_.u_i[j] = _.j_i[j] / _.rho_i;
			if (not inInitializationProcedure) {
				velocities.at(j)(i) = _.u_i[j];
			} else {
				_.u_i[j] = velocities.at(j)(i);
			}
		}
		densities(i) = _.rho_i;

		// apply collision operator
		collideOne<D, Q>(_);

		// copy back
		for (size_t j = 0; j < Q; j++) {
			F[j](i) = _.f_i[j];
		}

		// write omega2
		(*EntropicStabilized::OMEGA2)(i) = _.omega2;
		_.omega_mean += _.omega2;

	}
	//cout << _.omega_mean / locally_owned_dofs.size() << endl;

}
template void EntropicStabilized::collide<2, 9>(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure = false) const;
template void EntropicStabilized::collide<3, 19>(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure = false) const;
template void EntropicStabilized::collide<3, 27>(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure = false) const;


// ======================================================================================================

template<size_t D, size_t Q>
void EntropicStabilized::collideOne(EStCollisionData<D, Q>& _) const {
	/// PRECONDITIONS
	/// rho_i is set properly
	/// u_i is set properly
	/// tau is set
	/// f_post_i is the pre-collision distribution

	// BGK collision + regularization
	for (size_t j = 0; j < Q; j++) {
		_.f_post_i[j] = _.f_i[j];
	}
	//cout << _.f_post_i[0] << " " << _.f_i[0] << " " << _.rho_i << " " << _.u_i[0] << " " << _.u_i[1] << " " << _.tau << endl;
	collideBGK<D, Q>(_.f_post_i, _.rho_i, _.u_i, _.tau, _.stencil, _.f_eq_i);
	applyStabilizer<Q>(_.f_post_i, _.f_post_reg_i);
	//cout << "reg" << _.f_post_i[0] << _.f_post_reg_i[0] << endl;

	// calculate Kullback-Leibler Divergences (if f_i > 0), to determine omega2
	bool is_negative = false;
	for (size_t j = 0; j < Q; j++) {
		if (_.f_post_reg_i[j] < 1e-11)  {
			is_negative = true;
			break;
		}
		if (_.f_eq_i[j] < 1e-11)  {
			is_negative = true;
			break;
		}
		if (_.f_i[j] < 1e-11) {
			is_negative = true;
			break;
		}

	}
	if (is_negative) {
		//cout << "negative" << endl;
		_.omega2 = 1; // relaxation to stabilized post-collision state, without further manipulation of the mirror state
	} else {
		_.kld_pre = kullbackLeiblerDivergence<Q>(_.f_i, _.f_eq_i, _.rho_i);

		_.kld_post = kullbackLeiblerDivergence<Q>(_.f_post_reg_i, _.f_eq_i,
				_.rho_i);
		// alter relaxation so that the (stabilized) mirror distribution
		// has similar KLD as the pre-collision distribution,
		// i.e., the information from the high-order moments is incorporated
		// into the shear moments
		_.omega2 = _.kld_pre * abs(1.0 - 1. / _.tau) / _.kld_post;
		if (abs(_.kld_pre) < 1e-10){
			_.omega2 = 1;
		}
		//cout << _.kld_pre << " " <<  _.kld_post << " " << abs(1.0 - 1. / _.tau) << " " << _.omega2 << endl;
	}

	// constrain omega2. minimum: equilibrium distribution (0), maximum: full (over-)relaxation (2tau)
	// f^pc = f^eq + omega2 * (f^reg - f^eq)
	if (_.omega2 < 0.0){
		_.omega2 = 0.0;
	}
	else if (_.omega2 > 2 * _.tau){
		_.omega2 = 2 * _.tau;
	}

	// perform collision
	for (size_t j = 0; j < Q; j++) {
		_.f_i[j] = _.f_eq_i[j] + _.omega2 * (_.f_post_reg_i[j] - _.f_eq_i[j]);
	}
	//cout << "res" << _.f_i[0] << " " << _.f_eq_i[0] << " " << _.f_post_reg_i[0] << endl;
}
template void EntropicStabilized::collideOne<2,9>(EStCollisionData<2, 9>& _) const;
template void EntropicStabilized::collideOne<3,19>(EStCollisionData<3, 19>& _) const;
template void EntropicStabilized::collideOne<3,27>(EStCollisionData<3, 27>& _) const;

// ======================================================================================================

template<size_t Q>
double kullbackLeiblerDivergence(const array<double,Q>& f, const array<double,Q>& f_reg,
		double rho) {
	for (size_t i = 0; i < Q; i++) {
		assert(f[i] > 0.0);
		assert(f_reg[i] > 0.0);
	}
	double result = 0.0;
	for (size_t i = 0; i < Q; i++) {
		result += f[i] / rho * log(f[i] / f_reg[i]);
	}
	return result;
}
template double kullbackLeiblerDivergence<9>(const array<double,9>& f,
		const array<double,9>& f_reg, double rho);
template double kullbackLeiblerDivergence<19>(const array<double,19>& f,
		const array<double,19>& f_reg, double rho);
template double kullbackLeiblerDivergence<27>(const array<double,27>& f,
		const array<double,27>& f_reg, double rho);

// ======================================================================================================

template<size_t D, size_t Q>
void collideBGK(array<double,Q>& f, double rho, const array<double,D>& u, double tau,
		const Stencil& st, array<double,Q>& feq) {

	const vector<numeric_vector>& e = st.getDirections();
	const vector<double>& w = st.getWeights();
	double cs2 = st.getSpeedOfSoundSquare();

// calculate equilibrium distribution (feq)
	double uu_term = 0.0;
	for (size_t j = 0; j < D; j++) {
		uu_term += -(u[j] * u[j]) / (2.0 * cs2);
	}

	for (size_t i = 0; i < Q; i++) {
		double prefactor = w[i] * rho;

		double ue_term = 0.0;
		for (size_t j = 0; j < D; j++) {
			ue_term += (u[j] * e[i][j]) / st.getSpeedOfSoundSquare();
		}
		feq[i] = prefactor * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
	}
// update distribution
	for (size_t i = 0; i < Q; i++) {
		f[i] = f[i] - 1. / tau * (f[i] - feq[i]);
	}
}

template void collideBGK<2, 9>(array<double,9>& f, double rho, const array<double,2>& u, double tau,
		const Stencil& st, array<double,9>& feq);
template void collideBGK<3, 19>(array<double,19>& f, double rho, const array<double,3>& u,
		double tau, const Stencil& st, array<double,19>& feq);
template void collideBGK<3, 27>(array<double,27>& f, double rho, const array<double,3>& u,
		double tau, const Stencil& st, array<double,27>& feq);

// ======================================================================================================

template double density<9>(const array<double,9>& f);
template double density<19>(const array<double,19>& f);
template double density<27>(const array<double,27>& f);

template void momentum<2, 9>(const array<double,9>& f, array<double,2>& u,
		const vector<numeric_vector>& e);
template void momentum<3, 19>(const array<double,19>& f, array<double,3>& u,
		const vector<numeric_vector>& e);
template void momentum<3, 27>(const array<double,27>& f, array<double,3>& u,
		const vector<numeric_vector>& e);

// ======================================================================================================

template struct EStCollisionData<2, 9> ;
template struct EStCollisionData<3, 19> ;
template struct EStCollisionData<3, 27> ;

// ======================================================================================================


}
/* namespace natrium */

