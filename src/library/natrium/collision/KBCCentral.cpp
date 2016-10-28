/*
 * KBCCentral.cpp
 *
 *  Created on: 29.03.2016
 *      Author: Dominik Wilde
 */

#include "KBCCentral.h"
#define EVALUATE_GAMMA // if defined, an evaluation over time of the stabilizer gamma will be carried out

namespace natrium {

KBCCentral::KBCCentral(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil), m_D(setMRTWeights()), m_S(
				setRelaxationRates()) {

}

KBCCentral::~KBCCentral() {
}

void KBCCentral::collideAll(DistributionFunctions& f,
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
		throw CollisionException("KBC_Central only implemented for D2Q9/D3Q19");
	}
}

void KBCCentral::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	size_t Q = getQ();

	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare() / scaling / scaling;
	vector<double> meq(Q); 	// moment equilibrium distribution functions
	vector<double> m(Q);   	// moment distribution

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		// transform the velocity space into moment space
		m.at(0) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i)
				+ f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
		m.at(1) = -4 * f.at(0)(i) - f.at(1)(i) - f.at(2)(i) - f.at(3)(i)
				- f.at(4)(i)
				+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
//			m.at(2) = 4 * f.at(0)(i)
//					- 2 * (f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i))
		+f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
		m.at(3) = f.at(1)(i) - f.at(3)(i) + f.at(5)(i) - f.at(6)(i) - f.at(7)(i)
				+ f.at(8)(i);
//			m.at(4) = -2 * (f.at(1)(i) - f.at(3)(i)) + f.at(5)(i) - f.at(6)(i)
//					- f.at(7)(i) + f.at(8)(i);
		m.at(5) = f.at(2)(i) - f.at(4)(i) + f.at(5)(i) + f.at(6)(i) - f.at(7)(i)
				- f.at(8)(i);
//			m.at(6) = -2 * (f.at(2)(i) - f.at(4)(i)) + f.at(5)(i) + f.at(6)(i)
//					- f.at(7)(i) - f.at(8)(i);
		m.at(7) = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
		m.at(8) = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);

		// calculate the moment equilibrium distribution function
		double rho = m.at(0);
		double jx = m.at(3);
		double jy = m.at(5);

		// calculate density
		densities(i) = rho;

		if (not inInitializationProcedure) {

			velocities.at(0)(i) = scaling / densities(i) * jx;
			velocities.at(1)(i) = scaling / densities(i) * jy;
		}

		meq.at(0) = rho;
		meq.at(1) = -2 * rho + 3 / rho * (jx * jx + jy * jy);
		meq.at(2) = 0;
		meq.at(3) = jx;
		meq.at(4) = 0;
		meq.at(5) = jy;
		meq.at(6) = 0;
		meq.at(7) = jx * jx - jy * jy;
		meq.at(8) = jx * jy;

		//relax and rescale the moments

		m.at(1) = m.at(1)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(1) - meq.at(1));
		m.at(7) = m.at(7)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(7) - meq.at(7));
		m.at(8) = m.at(8)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(8) - meq.at(8));

//cout << getPrefactor() << " " << -1./(getRelaxationParameter()+0.5) << " " << getTime() << endl;

		m.at(2) = -rho - m.at(1);
		m.at(4) = -jx;
		m.at(6) = -jy;

		for (size_t j = 0; j < Q; j++) {
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
		f.at(5)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3) + m.at(4)
				+ m.at(5) + m.at(6) + m.at(8);
		f.at(6)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3) - m.at(4)
				+ m.at(5) + m.at(6) - m.at(8);
		f.at(7)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3) - m.at(4)
				- m.at(5) - m.at(6) + m.at(8);
		f.at(8)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3) + m.at(4)
				- m.at(5) - m.at(6) - m.at(8);

	}

//#ifdef EVALUATE_GAMMA
//	writeDeviation(gamma.getAverage(), gamma.getDeviation(),
//			entropy.getAverage(), entropy.getDeviation());
//#endif

}

void KBCCentral::collideAllD3Q19(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	double tm[19][19] = { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, { -30.0, -11.0,
			-11.0, -11.0, -11.0, -11.0, -11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
			8.0, 8.0, 8.0, 8.0, 8.0, 8.0 }, { 12.0, -4.0, -4.0, -4.0, -4.0,
			-4.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0 }, { 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0,
			1.0, -1.0, 1.0, -1.0, 0, 0, 0, 0 }, { 0.0, -4.0, 4.0, 0.0, 0.0, 0.0,
			0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0, 0, 0, 0 }, {
			0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0, 0, 0,
			0, 1.0, -1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0,
			1.0, 1.0, -1.0, -1.0, 0, 0, 0, 0, 1.0, -1.0, 1.0, -1.0 }, { 0.0,
			0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0,
			-1.0, 1.0, 1.0, -1.0, -1.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0,
			0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0 }, {
			0.0, 2.0, 2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0, 1.0, -2.0, -2.0, -2.0, -2.0 }, { 0.0, -4.0, -4.0, 2.0, 2.0,
			2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, -2.0,
			-2.0 }, { 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
			-1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, -2.0,
			-2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0,
			0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,
			-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			1.0, 0 - 1.0, -1.0, 1.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0,
			1.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			-1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0 }, {
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
			-1.0, -1.0, -1.0, -1.0, 1.0, 1.0 } };

	double invm[19][19]= {{0.0526316, -0.0125313, 0.047619, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0, 0, 0, 0, -0},
	{0.0526316, -0.00459482, -0.015873, 0.1, -0.1, -4.44089e-18, -1.11022e-18, 7.39557e-34, 1.84889e-34, 0.0555556, -0.0555556, -1.85037e-18, -9.25186e-19, 2.22045e-17, -0, -2.15704e-33, 1.11022e-17, 0, -0},
	{0.0526316, -0.00459482, -0.015873, -0.1, 0.1, -1.11022e-18, -2.77556e-19, 4.93038e-34, 1.2326e-34, 0.0555556, -0.0555556, 1.85037e-18, 9.25186e-19, 5.55112e-18, -0, -1.54074e-34, 2.77556e-18, -6.93889e-18, -0},
	{0.0526316, -0.00459482, -0.015873, -3.33067e-18, -8.32667e-19, 0.1, -0.1, -2.77556e-18, -6.93889e-19, -0.0277778, 0.0277778, 0.0833333, -0.0833333, -3.08149e-34, 6.93889e-18, 3.85186e-34, -6.93889e-18, 1.04083e-17, 3.46945e-18},
	{0.0526316, -0.00459482, -0.015873, 5.55112e-19, 1.38778e-19, -0.1, 0.1, -6.93889e-18, -8.67362e-18, -0.0277778, 0.0277778, 0.0833333, -0.0833333, 4.62223e-34, -1.04083e-17, 6.93889e-18, 1.04083e-17, -8.67362e-18, 1.73472e-18},
	{0.0526316, -0.00459482, -0.015873, 0, 0, 0, 0, 0.1, -0.1, -0.0277778, 0.0277778, -0.0833333, 0.0833333, 0, -0, 0, 0, 0, -6.93889e-18},
	{0.0526316, -0.00459482, -0.015873, 0, 0, 0, 0, -0.1, 0.1, -0.0277778, 0.0277778, -0.0833333, 0.0833333, 0, -0, 0, 0, 0, -0},
	{0.0526316, 0.00334169, 0.00396825, 0.1, 0.025, 0.1, 0.025, 2.22045e-17, 5.55112e-18, 0.0277778, 0.0138889, 0.0833333, 0.0416667, 0.25, 1.38778e-17, -6.93889e-18, 0.125, -0.125, -6.93889e-18},
	{0.0526316, 0.00334169, 0.00396825, -0.1, -0.025, 0.1, 0.025, 5.55112e-18, 1.38778e-18, 0.0277778, 0.0138889, 0.0833333, 0.0416667, -0.25, -0, 6.93889e-18, -0.125, -0.125, -0},
	{0.0526316, 0.00334169, 0.00396825, 0.1, 0.025, -0.1, -0.025, -2.77556e-18, -6.93889e-19, 0.0277778, 0.0138889, 0.0833333, 0.0416667, -0.25, 6.93889e-18, -6.93889e-18, 0.125, 0.125, -3.46945e-18},
	{0.0526316, 0.00334169, 0.00396825, -0.1, -0.025, -0.1, -0.025, 2.77556e-18, 6.93889e-19, 0.0277778, 0.0138889, 0.0833333, 0.0416667, 0.25, -6.93889e-18, 6.93889e-18, -0.125, 0.125, 3.46945e-18},
	{0.0526316, 0.00334169, 0.00396825, 0.1, 0.025, -4.44089e-18, -1.11022e-18, 0.1, 0.025, 0.0277778, 0.0138889, -0.0833333, -0.0416667, -5.55112e-18, -0, 0.25, -0.125, 0, 0.125},
	{0.0526316, 0.00334169, 0.00396825, -0.1, -0.025, -1.11022e-18, -2.77556e-19, 0.1, 0.025, 0.0277778, 0.0138889, -0.0833333, -0.0416667, 5.55112e-18, 1.38778e-17, -0.25, 0.125, -6.93889e-18, 0.125},
	{0.0526316, 0.00334169, 0.00396825, 0.1, 0.025, 1.11022e-18, 2.77556e-19, -0.1, -0.025, 0.0277778, 0.0138889, -0.0833333, -0.0416667, -5.55112e-18, -0, -0.25, -0.125, 6.93889e-18, -0.125},
	{0.0526316, 0.00334169, 0.00396825, -0.1, -0.025, 1.66533e-18, 4.16334e-19, -0.1, -0.025, 0.0277778, 0.0138889, -0.0833333, -0.0416667, 5.55112e-18, -1.38778e-17, 0.25, 0.125, -3.46945e-18, -0.125},
	{0.0526316, 0.00334169, 0.00396825, -3.33067e-18, -8.32667e-19, 0.1, 0.025, 0.1, 0.025, -0.0555556, -0.0277778, 0, 0, -3.08149e-34, 0.25, 3.85186e-34, -6.93889e-18, 0.125, -0.125},
	{0.0526316, 0.00334169, 0.00396825, 3.33067e-18, 8.32667e-19, -0.1, -0.025, 0.1, 0.025, -0.0555556, -0.0277778, 0, 0, 3.08149e-34, -0.25, -3.85186e-34, 6.93889e-18, -0.125, -0.125},
	{0.0526316, 0.00334169, 0.00396825, -3.33067e-18, -8.32667e-19, 0.1, 0.025, -0.1, -0.025, -0.0555556, -0.0277778, 0, 0, -3.08149e-34, -0.25, 3.85186e-34, -6.93889e-18, 0.125, 0.125},
	{0.0526316, 0.00334169, 0.00396825, 3.33067e-18, 8.32667e-19, -0.1, -0.025, -0.1, -0.025, -0.0555556, -0.0277778, 0, 0, 3.08149e-34, 0.25, -3.85186e-34, 6.93889e-18, -0.125, 0.125}};

	size_t Q = getQ();

	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare() / scaling / scaling;
	vector<double> meq(Q); 	// moment equilibrium distribution functions
	vector<double> m(Q);   	// moment distribution

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		for (int p = 0; p < 19; p++) {
			m[p]=0;
			for (int q = 0; q < 19; q++) {
				m[p] += tm[p][q] * f.at(q)(i);
			}
		}

		double rho = m.at(0);
		double jx = m.at(3);
		double jy = m.at(5);
		double jz = m.at(7);


		densities(i) = rho;

		if (not inInitializationProcedure) {

			velocities.at(0)(i) = scaling / densities(i) * jx;
			velocities.at(1)(i) = scaling / densities(i) * jy;
			velocities.at(2)(i) = scaling / densities(i) * jz;
		}

		meq.at(1) = -11 * rho + 19. / rho * (jx * jx + jy * jy + jz * jz);
		meq.at(9) = 1. / rho * (2 * jx * jx - (jy * jy + jz * jz));
		meq.at(11) = 1. / rho * (jy * jy - jz * jz);
		meq.at(13) = 1. / rho * jx * jy;
		meq.at(14) = 1. / rho * jy * jz;
		meq.at(15) = 1. / rho * jx * jz;

		m.at(1) = m.at(1)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(1) - meq.at(1));
		m.at(9) = m.at(9)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(9) - meq.at(9));
		m.at(11) = m.at(11)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(11) - meq.at(11));
		m.at(13) = m.at(13)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(13) - meq.at(13));
		m.at(14) = m.at(14)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(14) - meq.at(14));
		m.at(15) = m.at(15)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(15) - meq.at(15));


		m.at(2) = -7. / 38 * rho - 11. / 38 * m.at(1);
		m.at(4) = -2. / 3. * jx;
		m.at(6) = -2. / 3. * jy;
		m.at(8) = -2. / 3. * jz;
		m.at(10) = -1. / 2. * m.at(9);
		m.at(12) = -1. / 2. * m.at(11);
		m.at(16) = 0;
		m.at(17) = 0;
		m.at(18) = 0;

		for (int p = 0; p < 19; p++) {
			f.at(p)(i)=0;
			for (int q = 0; q < 19; q++) {
				f.at(p)(i) += invm[p][q] * m.at(q);
			}
		}

	}
}
}
