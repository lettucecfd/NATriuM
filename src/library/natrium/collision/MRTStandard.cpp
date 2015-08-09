/*
 * MRTStandard.cpp
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#include "MRTStandard.h"

namespace natrium {

MRTStandard::MRTStandard(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil) {

}

MRTStandard::~MRTStandard() {
}

double MRTStandard::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

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

}

double MRTStandard::getMomentEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

	double x;
	double scaling = getStencil()->getScaling();
	switch (i) {
	case 0: {
		x = rho;
		break;
	}
	case 1: {
		x = rho * (-2 + 3 * (u(0) * u(0) + u(1) * u(1)));
		break;
	}
	case 2: {
		x = rho * (1 - 3 * (u(0) * u(0) + u(1) * u(1)));
		break;
	}
	case 3: {
		x = rho * u(0);
		break;
	}
	case 4: {
		x = -rho * u(0);
		break;
	}
	case 5: {
		x = rho * u(1);
		break;
	}
	case 6: {
		x = -rho * u(1);
		break;
	}
	case 7: {
		x = rho * (u(0) * u(0) - u(1) * u(1));
		break;
	}
	case 8: {
		x = rho * u(0) * u(1);
		break;
	}
	default:
		x = 0;
	}
	return x;
}

void MRTStandard::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		bool inInitializationProcedure = false) const {
	size_t n_dofs = f.at(0).size();
	size_t Q = 9;
	size_t D = 2;
	double scaling = getStencil()->getScaling();
	vector<double> meq(Q); // moment equilibrium distribution functions
	vector<double> m(Q);

	for (size_t i = 0; i < n_dofs; i++) {

		velocities.at(0)(i) = scaling / densities(i)
				* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
						- f.at(6)(i) - f.at(7)(i));
		velocities.at(1)(i) = scaling / densities(i)
				* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
						- f.at(7)(i) - f.at(8)(i));

			// calculate density
			densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
			+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
			+ f.at(8)(i);

			m.at(0) = f.at(i)(0) + f.at(i)(1) + f.at(i)(2) + f.at(i)(3) + f.at(i)(4)
			+ f.at(i)(5) + f.at(i)(6) + f.at(i)(7) + f.at(i)(8);
			m.at(1) = -4 * f.at(i)(0) - f.at(i)(1) - f.at(i)(2) - f.at(i)(3)
			- f.at(i)(4)
			+ 2 * (f.at(i)(5) + f.at(i)(6) + f.at(i)(7) + f.at(i)(8));
			m.at(2) = 4 * f.at(i)(0)
			- 2 * (f.at(i)(1) + f.at(i)(2) + f.at(i)(3) + f.at(i)(4))
			+ f.at(i)(5) + f.at(i)(6) + f.at(i)(7) + f.at(i)(8);
			m.at(3) = f.at(i)(1) - f.at(i)(3) + f.at(i)(5) - f.at(i)(6) - f.at(i)(7)
			+ f.at(i)(8);
			m.at(4) = -2 * (f.at(i)(1) - f.at(i)(3)) + f.at(i)(5) - f.at(i)(6)
			- f.at(i)(7) + f.at(i)(8);
			m.at(5) = f.at(i)(2) - f.at(i)(4) + f.at(i)(5) + f.at(i)(6) - f.at(i)(7)
			- f.at(i)(8);
			m.at(6) = -2 * (f.at(i)(2) - f.at(i)(4)) + f.at(i)(5) + f.at(i)(6)
			- f.at(i)(7) - f.at(i)(8);
			m.at(7) = f.at(i)(1) - f.at(i)(2) + f.at(i)(3) - f.at(i)(4);
			m.at(8) = f.at(i)(5) - f.at(i)(6) + f.at(i)(7) - f.at(i)(8);

			for (size_t j = 0; j < Q; j++) {
				numeric_vector u(2);
				u(0) = velocities.at(0)(i);
				u(1) = velocities.at(0)(j);
				meq.at(j) = getMomentEquilibriumDistribution(j, u, densities(i));
			}

			f.at(i)(0) = m.at(0) - 4 * (m.at(1) - m.at(2));
			f.at(i)(1) = m.at(0) - m.at(1) - 2 * (m.at(2) + m.at(4)) + m.at(3) + m.at(7);
			f.at(i)(2) = m.at(0) - m.at(1) - 2 * (m.at(2) + m.at(6)) + m.at(5) - m.at(7);
			f.at(i)(3) = m.at(0) - m.at(1) - 2 * (m.at(2) - m.at(4)) - m.at(3) + m.at(7);
			f.at(i)(4) = m.at(0) - m.at(1) - 2 * (m.at(2) - m.at(6)) - m.at(5) - m.at(7);
			f.at(i)(5) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3) + m.at(4) + m.at(5) + m.at(6)
			+ m.at(8);
			f.at(i)(6) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3) - m.at(4) + m.at(5) + m.at(6)
			- m.at(8);
			f.at(i)(7) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3) - m.at(4) - m.at(5) - m.at(6)
			+ m.at(8);
			f.at(i)(8) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3) + m.at(4) - m.at(5) - m.at(6)
			- m.at(8);

		}
	}

}
// namespace natrium

