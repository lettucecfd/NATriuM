/*
 * KBCStandard.cpp
 *
 *  Created on: 17.11.2015
 *      Author: dominik
 */

#include "KBCStandard.h"

namespace natrium {

KBCStandard::KBCStandard(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil) {

}

KBCStandard::~KBCStandard() {

}

void KBCStandard::collideAll(DistributionFunctions& f,
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
		throw CollisionException("KBC only implemented for D2Q9");
		// Inefficient collision
		//BGK::collideAll(f, densities, velocities, locally_owned_dofs,
		//		inInitializationProcedure);
	}
}

void KBCStandard::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	size_t Q = getQ();

	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare() / (scaling * scaling);

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// calculate density
		densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
				+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
				+ f.at(8)(i);

		// density
		double rho = densities(i);

		if (not inInitializationProcedure) {

			velocities.at(0)(i) = scaling / rho
					* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
							- f.at(6)(i) - f.at(7)(i));
			velocities.at(1)(i) = scaling / rho
					* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
							- f.at(7)(i) - f.at(8)(i));
		}

		double ux = velocities.at(0)(i) / scaling;
		double uy = velocities.at(1)(i) / scaling;
		double scalar_product = ux * ux + uy * uy;

		// moment representation of the populations
		double T = 0, N = 0, Pi_xy = 0, Q_xyy = 0, Q_xxy = 0, Q_yxx = 0, A = 0;

		// calculate the trace of the pressure tensor at unit density (T)
		T = f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i)
				+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
		T = T / rho;

		// calculate the normal stress difference at unit density (N)
		N = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
		N = N / rho;

		// calculate the off-diagonal component of the pressure tensor at unit density
		Pi_xy = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);
		Pi_xy = Pi_xy / rho;

		// calculate moments of third order
		Q_xyy = f.at(5)(i) - f.at(6)(i) - f.at(7)(i) + f.at(8)(i);
		Q_xyy = Q_xyy / rho;

		Q_yxx = f.at(5)(i) + f.at(6)(i) - f.at(7)(i) - f.at(8)(i);
		Q_yxx = Q_yxx / rho;

		// calculate moments of forth order
		A = f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
		A = A / rho;

		// kinematic part vector k
		vector<double> k(Q);

		k.at(0) = rho;
		k.at(1) = 0.5 * rho * ux;
		k.at(3) = 0.5 * rho * -ux;
		k.at(2) = 0.5 * rho * uy;
		k.at(4) = 0.5 * rho * -uy;
		k.at(5) = 0;
		k.at(6) = 0;
		k.at(7) = 0;
		k.at(8) = 0;

		// shear part vector s
		vector<double> s(Q);

		s.at(0) = -rho * T;
		s.at(1) = 0.5 * rho * 0.5 * (T + N);
		s.at(2) = 0.5 * rho * 0.5 * (T - N);
		s.at(3) = 0.5 * rho * 0.5 * (T + N);
		s.at(4) = 0.5 * rho * 0.5 * (T - N);
		s.at(5) = 0.25 * rho * Pi_xy;
		s.at(6) = 0.25 * rho * -Pi_xy;
		s.at(7) = 0.25 * rho * Pi_xy;
		s.at(8) = 0.25 * rho * -Pi_xy;

		// higher order part vector h
		vector<double> h(Q);

		h.at(0) = rho * A;
		h.at(1) = -0.5 * rho * (1 * Q_xyy + A);
		h.at(3) = -0.5 * rho * (-1 * Q_xyy + A);
		h.at(2) = -0.5 * rho * (1 * Q_yxx + A);
		h.at(4) = -0.5 * rho * (-1 * Q_yxx + A);
		h.at(5) = 0.25 * rho * (1 * Q_xyy + 1 * Q_yxx + A);
		h.at(6) = 0.25 * rho * (-1 * Q_xyy + 1 * Q_yxx + A);
		h.at(7) = 0.25 * rho * (-1 * Q_xyy - 1 * Q_yxx + A);
		h.at(8) = 0.25 * rho * (1 * Q_xyy + -1 * Q_yxx + A);

		// equilibrium vectors for shear vector
		vector<double> seq(Q);
		// equilibrium vectors for higher order part
		vector<double> heq(Q);

		// equilibrium distribution of the population f
		vector<double> feq(Q, 0.0);
		double uSquareTerm;
		double mixedTerm;
		double weighting;
		double prefactor = scaling / getStencil()->getSpeedOfSoundSquare();

		uSquareTerm = -scalar_product
				/ (2 * getStencil()->getSpeedOfSoundSquare());
		// direction 0
		weighting = 4. / 9. * rho;
		feq.at(0) = weighting * (1 + uSquareTerm);
		// directions 1-4
		weighting = 1. / 9. * rho;
		mixedTerm = prefactor * (velocities.at(0)(i));
		feq.at(1) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(3) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (velocities.at(1)(i));
		feq.at(2) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(4) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		// directions 5-8
		weighting = 1. / 36. * rho;
		mixedTerm = prefactor * (velocities.at(0)(i) + velocities.at(1)(i));
		feq.at(5) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(7) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (-velocities.at(0)(i) + velocities.at(1)(i));
		feq.at(6) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(8) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);

		//Calculation of the moments for the equilibrium distribution function

		//calculate the trace of the pressure tensor at unit density (T)
		T = feq.at(1) + feq.at(2) + feq.at(3) + feq.at(4)
				+ 2 * (feq.at(5) + feq.at(6) + feq.at(7) + feq.at(8));
		T = T / rho;

		//calculate the normal stress difference at unit density (N)
		N = feq.at(1) - feq.at(2) + feq.at(3) - feq.at(4);
		N = N / rho;

		//calculate the off-diagonal component of the pressure tensor at unit density
		Pi_xy = feq.at(5) - feq.at(6) + feq.at(7) - feq.at(8);
		Pi_xy = Pi_xy / rho;

		//calculate moments of third order
		Q_xyy = feq.at(5) - feq.at(6) - feq.at(7) + feq.at(8);
		Q_xyy = Q_xyy / rho;

		Q_yxx = feq.at(5) + feq.at(6) - feq.at(7) - feq.at(8);
		Q_yxx = Q_yxx / rho;

		// calculate moments of forth order
		A = feq.at(5) + feq.at(6) + feq.at(7) + feq.at(8);
		A = A / rho;

		// calculate shear equilibrium
		seq.at(0) = -rho * T;
		seq.at(1) = 0.5 * rho * 0.5 * (T + N);
		seq.at(3) = 0.5 * rho * 0.5 * (T + N);
		seq.at(2) = 0.5 * rho * 0.5 * (T - N);
		seq.at(4) = 0.5 * rho * 0.5 * (T - N);
		seq.at(5) = 0.25 * rho * Pi_xy;
		seq.at(6) = 0.25 * rho * -Pi_xy;
		seq.at(7) = 0.25 * rho * Pi_xy;
		seq.at(8) = 0.25 * rho * -Pi_xy;

		// calculate higher order equilibrium
		heq.at(0) = rho * A;
		heq.at(1) = -0.5 * rho * (1 * Q_xyy + A);
		heq.at(3) = -0.5 * rho * (-1 * Q_xyy + A);
		heq.at(2) = -0.5 * rho * (1 * Q_yxx + A);
		heq.at(4) = -0.5 * rho * (-1 * Q_yxx + A);
		heq.at(5) = 0.25 * rho * (1 * Q_xyy + 1 * Q_yxx + A);
		heq.at(6) = 0.25 * rho * (-1 * Q_xyy + 1 * Q_yxx + A);
		heq.at(7) = 0.25 * rho * (-1 * Q_xyy - 1 * Q_yxx + A);
		heq.at(8) = 0.25 * rho * (1 * Q_xyy + -1 * Q_yxx + A);

		//deviation of the shear parts
		vector<double> delta_s(Q);
		//deviation of the higher order parts
		vector<double> delta_h(Q);

		//entropic scalar product of s
		vector<double> delta_seq(Q);
		//entropic scalar product of h
		vector<double> delta_heq(Q);

		// needed expressions for the calculation of gamma
		delta_s.at(0) = s.at(0) - seq.at(0);
		delta_s.at(1) = s.at(1) - seq.at(1);
		delta_s.at(2) = s.at(2) - seq.at(2);
		delta_s.at(3) = s.at(3) - seq.at(3);
		delta_s.at(4) = s.at(4) - seq.at(4);
		delta_s.at(5) = s.at(5) - seq.at(5);
		delta_s.at(6) = s.at(6) - seq.at(6);
		delta_s.at(7) = s.at(7) - seq.at(7);
		delta_s.at(8) = s.at(8) - seq.at(8);

		delta_h.at(0) = h.at(0) - heq.at(0);
		delta_h.at(1) = h.at(1) - heq.at(1);
		delta_h.at(2) = h.at(2) - heq.at(2);
		delta_h.at(3) = h.at(3) - heq.at(3);
		delta_h.at(4) = h.at(4) - heq.at(4);
		delta_h.at(5) = h.at(5) - heq.at(5);
		delta_h.at(6) = h.at(6) - heq.at(6);
		delta_h.at(7) = h.at(7) - heq.at(7);
		delta_h.at(8) = h.at(8) - heq.at(8);

		delta_seq.at(0) = delta_s.at(0) * delta_h.at(0) / feq.at(0);
		delta_seq.at(1) = delta_s.at(1) * delta_h.at(1) / feq.at(1);
		delta_seq.at(2) = delta_s.at(2) * delta_h.at(2) / feq.at(2);
		delta_seq.at(3) = delta_s.at(3) * delta_h.at(3) / feq.at(3);
		delta_seq.at(4) = delta_s.at(4) * delta_h.at(4) / feq.at(4);
		delta_seq.at(5) = delta_s.at(5) * delta_h.at(5) / feq.at(5);
		delta_seq.at(6) = delta_s.at(6) * delta_h.at(6) / feq.at(6);
		delta_seq.at(7) = delta_s.at(7) * delta_h.at(7) / feq.at(7);
		delta_seq.at(8) = delta_s.at(8) * delta_h.at(8) / feq.at(8);

		delta_heq.at(0) = delta_h.at(0) * delta_h.at(0) / feq.at(0);
		delta_heq.at(1) = delta_h.at(1) * delta_h.at(1) / feq.at(1);
		delta_heq.at(2) = delta_h.at(2) * delta_h.at(2) / feq.at(2);
		delta_heq.at(3) = delta_h.at(3) * delta_h.at(3) / feq.at(3);
		delta_heq.at(4) = delta_h.at(4) * delta_h.at(4) / feq.at(4);
		delta_heq.at(5) = delta_h.at(5) * delta_h.at(5) / feq.at(5);
		delta_heq.at(6) = delta_h.at(6) * delta_h.at(6) / feq.at(6);
		delta_heq.at(7) = delta_h.at(7) * delta_h.at(7) / feq.at(7);
		delta_heq.at(8) = delta_h.at(8) * delta_h.at(8) / feq.at(8);

		// sum of all entropic scalar products of s
		double sum_s = 0;
		sum_s = delta_seq.at(0) + delta_seq.at(1) + delta_seq.at(2)
				+ delta_seq.at(3) + delta_seq.at(4) + delta_seq.at(5)
				+ delta_seq.at(6) + delta_seq.at(7) + delta_seq.at(8);

		// sum of all entropic scalar products of h
		double sum_h = 0;
		sum_h = delta_heq.at(0) + delta_heq.at(1) + delta_heq.at(2)
				+ delta_heq.at(3) + delta_heq.at(4) + delta_heq.at(5)
				+ delta_heq.at(6) + delta_heq.at(7) + delta_heq.at(8);

		// relaxation parameter (2*beta = omega of BGK)
		double beta = 1. / (getRelaxationParameter() + 0.5) / 2;

		// stabilizer of KBC model
		double gamma = 1. / beta - (2 - 1. / beta) * sum_s / sum_h;

		// if the sum_h expression is too small, BGK shall be performed (gamma = 2)
		if (sum_h < 1e-20) {
			gamma = 2;
		}

		// calculate new f
		f.at(0)(i) = f.at(0)(i)
				- beta * (2 * delta_s.at(0) + gamma * delta_h.at(0));
		f.at(1)(i) = f.at(1)(i)
				- beta * (2 * delta_s.at(1) + gamma * delta_h.at(1));
		f.at(2)(i) = f.at(2)(i)
				- beta * (2 * delta_s.at(2) + gamma * delta_h.at(2));
		f.at(3)(i) = f.at(3)(i)
				- beta * (2 * delta_s.at(3) + gamma * delta_h.at(3));
		f.at(4)(i) = f.at(4)(i)
				- beta * (2 * delta_s.at(4) + gamma * delta_h.at(4));
		f.at(5)(i) = f.at(5)(i)
				- beta * (2 * delta_s.at(5) + gamma * delta_h.at(5));
		f.at(6)(i) = f.at(6)(i)
				- beta * (2 * delta_s.at(6) + gamma * delta_h.at(6));
		f.at(7)(i) = f.at(7)(i)
				- beta * (2 * delta_s.at(7) + gamma * delta_h.at(7));
		f.at(8)(i) = f.at(8)(i)
				- beta * (2 * delta_s.at(8) + gamma * delta_h.at(8));

		if (gamma != 2) {

		}

	}

}
}
/* namespace natrium */
