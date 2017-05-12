/*
 * EntropicStabilized.h
 *
 *  Created on: 11.05.2017
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_COLLISION_ENTROPICSTABILIZED_H_
#define LIBRARY_NATRIUM_COLLISION_ENTROPICSTABILIZED_H_

#include "MRT.h"
#include <array>

using std::array;

namespace natrium {

template<size_t D, size_t Q>
struct EStCollisionData {
	const Stencil& stencil;
	const vector<numeric_vector>& e;
	double scaling;
	double tau;

	double omega2;
	double omega_mean;

	array<double,Q> f_i;
	array<double,Q> f_post_i;
	array<double,Q> f_post_reg_i;
	array<double,Q> f_eq_i;

	double kld_pre;
	double kld_post;
    double kld_post_reg;

	double rho_i;
	array<double,D> u_i;
	array<double,D> j_i;

	EStCollisionData(const Stencil& st, double tau) :
			stencil(st), e(st.getDirections()), scaling(st.getScaling()), tau(
					tau), omega2(0.0), omega_mean(0.0), kld_pre(0.0), kld_post(0.0), rho_i(0.0){
		assert(st.getQ() == Q);
		assert(st.getD() == D);

		/*double f_i[Q] = { };
		double f_post_i[Q] = { };
		double f_post_reg_i[Q] = { };
		double f_eq_i[Q] = { };
		double u_i[D] = { };
		double j_i[D] = { };*/

	}
};

class EntropicStabilized: public MRT {
public:
	EntropicStabilized(double relaxationParameter, double dt,
			const boost::shared_ptr<Stencil> stencil);
	virtual ~EntropicStabilized();

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	template<size_t D, size_t Q>
	void collide(DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	template<size_t D, size_t Q>
	void collideOne(EStCollisionData<D,Q>& _) const;

	static boost::shared_ptr<distributed_vector> OMEGA2;
};

template<size_t Q>
double kullbackLeiblerDivergence(const array<double,Q>& f, const array<double,Q>& f_reg,
		double rho);

template<size_t D, size_t Q>
void collideBGK(array<double,Q>& f, double rho, const array<double,D>& u, double tau,
		const Stencil& st, array<double,Q>& feq);

template<size_t Q>
double density(const array<double,Q>& f) {
	double result = 0.0;
	for (size_t i = 0; i < Q; i++) {
		result += f[i];
	}
	return result;
}

template<size_t D, size_t Q>
void momentum(const array<double,Q>& f, array<double,D>& j, const vector<numeric_vector>& e) {
	for (size_t i = 0; i < D; i++) {
		j[i] = 0.0;
	}
	for (size_t i = 0; i < D; i++) {
		for (size_t k = 0; k < Q; k++) {
			j[i] += f[k] * e[k](i);
		}
	}
}

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_COLLISION_ENTROPICSTABILIZED_H_ */
