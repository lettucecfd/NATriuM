/*
 * D2Q9PseudopotentialModel.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: kk
 */

#include "D2Q9PseudopotentialModel.h"

#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

/// constructor
D2Q9PseudopotentialModel::D2Q9PseudopotentialModel(double scaling,
		const double dt) :
		D2Q9Model(scaling), m_dt(dt) {

}

// destructor
D2Q9PseudopotentialModel::~D2Q9PseudopotentialModel() {
}

// getEquilibriumDistribution
double D2Q9PseudopotentialModel::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {
	assert(i < Q);
	assert(rho > 0);
	assert(i >= 0);
	assert(u.size() == D);
	assert(u(0) < 1000000000000000.);
	assert(u(1) < 1000000000000000.);

	double prefactor = getWeight(i) * rho;
	double uSquareTerm = -Math::scalar_product(u, u)
			/ (2 * m_speedOfSoundSquare);
	if (0 == i) {
		return prefactor * (1 + uSquareTerm);
	}
	double mixedTerm = Math::scalar_product(u, getDirection(i))
			/ m_speedOfSoundSquare;
	return prefactor * (1 + mixedTerm * (1 + 0.5 * (mixedTerm)) + uSquareTerm);
}

// getInteractionForce
void D2Q9PseudopotentialModel::getInteractionForce(
		const vector<double>& distributions, numeric_vector & interactionForce,
		const double rho) {

	const double G = -4.4;
	numeric_vector densityGradient(2);

	interactionForce(0) = 0.0;
	interactionForce(1) = 0.0;

	//TODO getDensityGradient ;

	interactionForce(0) = -G * (1. - exp(rho)) * densityGradient(0);
	interactionForce(1) = -G * (1. - exp(rho)) * densityGradient(1);
}
}
