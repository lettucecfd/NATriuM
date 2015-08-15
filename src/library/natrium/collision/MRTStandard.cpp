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
		MRT(relaxationParameter, dt, stencil), m_D(setMRTWeights()), m_S(
				setRelaxationRates()) {

}

MRTStandard::~MRTStandard() {
}

void MRTStandard::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		bool inInitializationProcedure = false) const {

	if (Stencil_D2Q9 != getStencil()->getStencilType()) {
		// No MRT collision for other than D2Q9
		throw CollisionException("MRT only implemented for D2Q9");
	} else {

		size_t n_dofs = f.at(0).size();
		size_t Q = getQ();
		double scaling = getStencil()->getScaling();
		vector<double> meq(Q); 	// moment equilibrium distribution functions
		vector<double> m(Q);   	// moment distribution

		for (size_t i = 0; i < n_dofs; i++) {

			if (not inInitializationProcedure) {

				velocities.at(0)(i) = scaling / densities(i)
						* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
								- f.at(6)(i) - f.at(7)(i));
				velocities.at(1)(i) = scaling / densities(i)
						* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
								- f.at(7)(i) - f.at(8)(i));
			}
			// calculate density
			densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);
			if (i == 1) {
			}

			// transform the velocity space into moment space
			m.at(0) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);
			m.at(1) = -4 * f.at(0)(i) - f.at(1)(i) - f.at(2)(i) - f.at(3)(i)
					- f.at(4)(i)
					+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
			m.at(2) = 4 * f.at(0)(i)
					- 2 * (f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i))
					+ f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
			m.at(3) = f.at(1)(i) - f.at(3)(i) + f.at(5)(i) - f.at(6)(i)
					- f.at(7)(i) + f.at(8)(i);
			m.at(4) = -2 * (f.at(1)(i) - f.at(3)(i)) + f.at(5)(i) - f.at(6)(i)
					- f.at(7)(i) + f.at(8)(i);
			m.at(5) = f.at(2)(i) - f.at(4)(i) + f.at(5)(i) + f.at(6)(i)
					- f.at(7)(i) - f.at(8)(i);
			m.at(6) = -2 * (f.at(2)(i) - f.at(4)(i)) + f.at(5)(i) + f.at(6)(i)
					- f.at(7)(i) - f.at(8)(i);
			m.at(7) = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
			m.at(8) = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);

			// calculate the moment equilibrium distribution function
			meq.at(0) = densities(i);
			meq.at(1) = densities(i)
					* (-2
							+ 3
									* (velocities.at(0)(i) * velocities.at(0)(i)
											+ velocities.at(1)(i)
													* velocities.at(1)(i)));
			meq.at(2) = densities(i)
					* (1
							- 3
									* (velocities.at(0)(i) * velocities.at(0)(i)
											+ velocities.at(1)(i)
													* velocities.at(1)(i)));
			meq.at(3) = densities(i) * velocities.at(0)(i);
			meq.at(4) = -densities(i) * velocities.at(0)(i);
			meq.at(5) = densities(i) * velocities.at(1)(i);
			meq.at(6) = -densities(i) * velocities.at(1)(i);
			meq.at(7) = densities(i)
					* (velocities.at(0)(i) * velocities.at(0)(i)
							- velocities.at(1)(i) * velocities.at(1)(i));
			meq.at(8) = densities(i) * velocities.at(0)(i)
					* velocities.at(1)(i);

			//relax and rescale the moments
			for (size_t j = 0; j < Q; j++) {
				m.at(j) = m.at(j) - m_S.at(j) * (m.at(j) - meq.at(j));
				m.at(j) /= m_D.at(j);
			}

			//transform the momentum space back into velocity space
			f.at(0)(i) = m.at(0) - 4 * (m.at(1) - m.at(2));
			f.at(1)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) + m.at(4)) + m.at(3)
					+ m.at(7);
			f.at(2)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) + m.at(6)) + m.at(5)
					- m.at(7);
			f.at(3)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) - m.at(4)) - m.at(3)
					+ m.at(7);
			f.at(4)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) - m.at(6)) - m.at(5)
					- m.at(7);
			f.at(5)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3)
					+ m.at(4) + m.at(5) + m.at(6) + m.at(8);
			f.at(6)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3)
					- m.at(4) + m.at(5) + m.at(6) - m.at(8);
			f.at(7)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3)
					- m.at(4) - m.at(5) - m.at(6) + m.at(8);
			f.at(8)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3)
					+ m.at(4) - m.at(5) - m.at(6) - m.at(8);

		}
	}
}   	// if - else
}
// namespace natrium

