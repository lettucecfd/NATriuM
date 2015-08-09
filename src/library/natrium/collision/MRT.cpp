/*
 * MRT.cpp
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#include "MRT.h"

namespace natrium {

MRT::MRT(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		CollisionModel(stencil), m_relaxationParameter(relaxationParameter), m_prefactor(
				-1. / (relaxationParameter + 0.5)), m_dt(dt) {

}

MRT::~MRT() {
	// TODO Auto-generated destructor stub
}
/*
void MRT::collideSinglePoint(vector<double>& distributions) const {

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
} */

} // namespace natrium

