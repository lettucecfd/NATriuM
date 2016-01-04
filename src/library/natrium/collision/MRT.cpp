/*
 * MRT.cpp
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#include "MRT.h"

namespace natrium {

MRT::MRT(double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil) :
		CollisionModel(stencil), m_relaxationParameter(relaxationParameter), m_prefactor(
				-1. / (relaxationParameter + 0.5)), m_dt(dt) {

}

MRT::~MRT() {
	// TODO Auto-generated destructor stub
}

double MRT::getEquilibriumDistribution(size_t i, const numeric_vector& u,
		const double rho) const {

	assert(i < getStencil()->getQ());
	assert(rho > 0);
	assert(i >= 0);
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

} // The BGK Equilibrium Distribution is used to initialize the system in the beginning

} // namespace natrium

