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
		double rho = densities(i);

		if (not inInitializationProcedure) {

			velocities.at(0)(i) = scaling / densities(i)
					* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
							- f.at(6)(i) - f.at(7)(i));
			velocities.at(1)(i) = scaling / densities(i)
					* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
							- f.at(7)(i) - f.at(8)(i));
		}

		double ux = velocities.at(0)(i) / scaling;
		double uy = velocities.at(1)(i) / scaling;
		double scalar_product = ux * ux + uy * uy;

		//moment representation of the populations
		double T = 0, N = 0, Pi_xy = 0, Q_xyy = 0, Q_yxx = 0, A = 0;

		//calculate the trace of the pressure tensor at unit density (T)
		T = f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i)
				+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
		T = T / densities(i);

		//calculate the normal stress difference at unit density (N)
		N = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
		N = N / densities(i);

		//calculate the off-diagonal component of the pressure tensor at unit density
		Pi_xy = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);
		Pi_xy = Pi_xy / densities(i);

		//calculate moments of third order
		Q_xyy = f.at(5)(i) - f.at(6)(i) - f.at(7)(i) + f.at(8)(i);
		Q_xyy = Q_yxx / densities(i);

		Q_yxx = f.at(5)(i) + f.at(6)(i) - f.at(7)(i) - f.at(8)(i);
		Q_yxx = Q_yxx / densities(i);

		//calculate moments of forth order
		A = f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
		A = A / densities(i);

		//transform from natural into central moments
		T = T - scalar_product;

		N = N - (ux * ux - uy * uy);

		Pi_xy = Pi_xy - ux * uy;

		Q_xyy = Q_xyy - 2 * uy * Pi_xy + 0.5 * ux * N - 0.5 * ux * T
				- ux * uy * uy;
		Q_yxx = Q_yxx - 2 * ux * Pi_xy - 0.5 * uy * N - 0.5 * uy * T
				- uy * ux * ux;

		A = A - 2 * (ux * Q_xyy + uy * Q_yxx) - 4 * ux * uy * Pi_xy
				- 0.5 * scalar_product * T + 0.5 * (ux * ux - uy * uy) * N
				- ux * ux * uy * uy;

		// kinematic part vector k
		vector<double> k(Q);

		k.at(0) = densities(i) * (1 - scalar_product);
		k.at(1) = 0.5 * densities(i) * (ux * ux + 1 * ux);
		k.at(3) = 0.5 * densities(i) * (ux * ux - 1 * ux);
		k.at(2) = 0.5 * densities(i) * (uy * uy + 1 * uy);
		k.at(4) = 0.5 * densities(i) * (uy * uy - 1 * uy);
		k.at(5) = 0.25 * densities(i) * (1 * 1) * ux * uy;
		k.at(6) = 0.25 * densities(i) * (1 * -1) * ux * uy;
		k.at(7) = 0.25 * densities(i) * (-1 * -1) * ux * uy;
		k.at(8) = 0.25 * densities(i) * (-1 * 1) * ux * uy;

		//shear part vector s
		vector<double> s(Q);

		s.at(0) = rho
				* ((4 * ux) * uy * Pi_xy - 0.5 * (ux * ux - uy * uy) * N
						+ 0.5 * (scalar_product - 2) * T);
		s.at(1) = 0.5 * rho
				* (0.5 * (1 + 1 * ux + ux * ux - uy * uy) * N
						- (2 * 1 * uy + 4 * ux * uy) * Pi_xy
						+ 0.5 * (1 - 1 * ux - scalar_product) * T);

		s.at(3) = 0.5 * rho
				* (0.5 * (1 - 1 * ux + ux * ux - uy * uy) * N
						- (2 * -1 * uy + 4 * ux * uy) * Pi_xy
						+ 0.5 * (1 + 1 * ux - scalar_product) * T);

		s.at(2) = 0.5 * rho
				* (0.5 * (-1 - 1 * uy + ux * ux - uy * uy) * N
						- (2 * 1 * ux + 4 * ux * uy) * Pi_xy
						+ 0.5 * (1 - 1 * uy - scalar_product) * T);

		s.at(4) = 0.5 * rho
				* (0.5 * (-1 + 1 * uy + ux * ux - uy * uy) * N
						- (2 * -1 * ux + 4 * ux * uy) * Pi_xy
						+ 0.5 * (1 + 1 * uy - scalar_product) * T);

		s.at(5) = 0.25 * rho
				* ((4 * ux * uy + (1) * (1) + 2 * 1 * uy + 2 * 1 * ux) * Pi_xy
						+ 0.5 * (-ux * ux + uy * uy - 1 * ux + 1 * uy) * N
						+ 0.5 * (scalar_product + 1 * ux + 1 * uy) * T);
		s.at(6) = 0.25 * rho
				* ((4 * ux * uy + (-1) * (1) + 2 * -1 * uy + 2 * 1 * ux) * Pi_xy
						+ 0.5 * (-ux * ux + uy * uy + 1 * ux + 1 * uy) * N
						+ 0.5 * (scalar_product + (-1) * ux + 1 * uy) * T);
		s.at(7) = 0.25 * rho
				* ((4 * ux * uy + (-1) * (-1) + 2 * -1 * uy + 2 * (-1) * ux)
						* Pi_xy
						+ 0.5 * (-ux * ux + uy * uy + 1 * ux - 1 * uy) * N
						+ 0.5 * (scalar_product + (-1) * ux - 1 * uy) * T);

		s.at(8) = 0.25 * rho
				* ((4 * ux * uy + (1) * (-1) + 2 * 1 * uy + 2 * (-1) * ux)
						* Pi_xy
						+ 0.5 * (-ux * ux + uy * uy - 1 * ux - 1 * uy) * N
						+ 0.5 * (scalar_product + (1) * ux - 1 * uy) * T);

		//higher order part vector h
		vector<double> h(Q);

		h.at(0) = rho * (2 * ux * Q_xyy + 2 * uy * Q_yxx + A);
		h.at(1) = 0.5 * rho * (-(1 + 2 * ux) * Q_xyy - 2 * uy * Q_yxx - A);
		h.at(3) = 0.5 * rho * (-(-1 + 2 * ux) * Q_xyy - 2 * uy * Q_yxx - A);
		h.at(2) = 0.5 * rho * (-(1 + 2 * uy) * Q_yxx - 2 * ux * Q_xyy - A);
		h.at(4) = 0.5 * rho * (-(-1 + 2 * uy) * Q_yxx - 2 * ux * Q_xyy - A);
		h.at(5) = 0.25 * rho
				* ((1 + 2 * ux) * Q_xyy + (1 + 2 * uy) * Q_yxx + A);
		h.at(6) = 0.25 * rho
				* ((-1 + 2 * ux) * Q_xyy + (1 + 2 * uy) * Q_yxx + A);
		h.at(7) = 0.25 * rho
				* ((-1 + 2 * ux) * Q_xyy + (-1 + 2 * uy) * Q_yxx + A);
		h.at(8) = 0.25 * rho
				* ((1 + 2 * ux) * Q_xyy + (-1 + 2 * uy) * Q_yxx + A);

		//equilibrium vectors for shear and higher order part vectors
		vector<double> seq(Q);
		vector<double> heq(Q);

		//equilibrium for T and A
		T = 2 * cs2;
		A = cs2 * cs2;

		//calculate shear equilibrium
		seq.at(0) = rho * (0.5 * (scalar_product - 2) * T);
		seq.at(1) = 0.5 * rho * 0.5 * (1 - ux - scalar_product) * T;
		seq.at(3) = 0.5 * rho * 0.5 * (1 + ux - scalar_product) * T;
		seq.at(2) = 0.5 * rho * 0.5 * (1 - uy - scalar_product) * T;
		seq.at(4) = 0.5 * rho * 0.5 * (1 + uy - scalar_product) * T;
		seq.at(5) = 0.25 * rho * 0.5 * (scalar_product + ux + uy) * T;
		seq.at(6) = 0.25 * rho * 0.5 * (scalar_product - ux + uy) * T;
		seq.at(7) = 0.25 * rho * 0.5 * (scalar_product - ux - uy) * T;
		seq.at(8) = 0.25 * rho * 0.5 * (scalar_product + ux - uy) * T;

		//calculate higher order equilibrium
		heq.at(0) = rho * A;
		heq.at(1) = -0.5 * rho * A;
		heq.at(3) = -0.5 * rho * A;
		heq.at(2) = -0.5 * rho * A;
		heq.at(4) = -0.5 * rho * A;
		heq.at(5) = 0.25 * rho * A;
		heq.at(6) = 0.25 * rho * A;
		heq.at(7) = 0.25 * rho * A;
		heq.at(8) = 0.25 * rho * A;

		// equilibrium distribution of the population f
		vector<double> feq(Q, 0.0);
		double uSquareTerm;
		double mixedTerm;
		double weighting;
		double prefactor = scaling / getStencil()->getSpeedOfSoundSquare();

		uSquareTerm = -scalar_product
				/ (2 * getStencil()->getSpeedOfSoundSquare());
		// direction 0
		weighting = 4. / 9. * densities(i);
		feq.at(0) = weighting * (1 + uSquareTerm);
		// directions 1-4
		weighting = 1. / 9. * densities(i);
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
		weighting = 1. / 36. * densities(i);
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

		vector<double> delta_s(Q);
		vector<double> delta_h(Q);

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

		delta_s.at(0) = delta_s.at(0) * delta_h.at(0) / feq.at(0);
		delta_s.at(1) = delta_s.at(1) * delta_h.at(1) / feq.at(1);
		delta_s.at(2) = delta_s.at(2) * delta_h.at(2) / feq.at(2);
		delta_s.at(3) = delta_s.at(3) * delta_h.at(3) / feq.at(3);
		delta_s.at(4) = delta_s.at(4) * delta_h.at(4) / feq.at(4);
		delta_s.at(5) = delta_s.at(5) * delta_h.at(5) / feq.at(5);
		delta_s.at(6) = delta_s.at(6) * delta_h.at(6) / feq.at(6);
		delta_s.at(7) = delta_s.at(7) * delta_h.at(7) / feq.at(7);
		delta_s.at(8) = delta_s.at(8) * delta_h.at(8) / feq.at(8);

		delta_h.at(0) = delta_h.at(0) * delta_h.at(0) / feq.at(0);
		delta_h.at(1) = delta_h.at(1) * delta_h.at(1) / feq.at(1);
		delta_h.at(2) = delta_h.at(2) * delta_h.at(2) / feq.at(2);
		delta_h.at(3) = delta_h.at(3) * delta_h.at(3) / feq.at(3);
		delta_h.at(4) = delta_h.at(4) * delta_h.at(4) / feq.at(4);
		delta_h.at(5) = delta_h.at(5) * delta_h.at(5) / feq.at(5);
		delta_h.at(6) = delta_h.at(6) * delta_h.at(6) / feq.at(6);
		delta_h.at(7) = delta_h.at(7) * delta_h.at(7) / feq.at(7);
		delta_h.at(8) = delta_h.at(8) * delta_h.at(8) / feq.at(8);

		double sum_h = delta_h.at(0) + delta_h.at(1) + delta_h.at(2)
				+ delta_h.at(3) + delta_h.at(4) + delta_h.at(5) + delta_h.at(6)
				+ delta_h.at(7) + delta_h.at(8);

		double sum_s = delta_s.at(0) + delta_s.at(1) + delta_s.at(2)
				+ delta_s.at(3) + delta_s.at(4) + delta_s.at(5) + delta_s.at(6)
				+ delta_s.at(7) + delta_s.at(8);

		//	double viscosity = getRelaxationParameter() * cs2;

		double beta = 1. / (getRelaxationParameter() + 0.5) / 2.; // getestet und fuer gut befunden

		double gamma = 1. / beta - (2 - 1. / beta) * sum_s / sum_h;

		if (sum_h < 1e-15) {
			//gamma = 1. / beta - (2 - 1. / beta);
			gamma = 2.0;
		}

		gamma = 2.0; //only for testing;
		for (int z = 0; z < 9; z++) {
		 //assert (fabs (k.at(z)+s.at(z) + h.at(z) - f.at(z)(i) ) < 1e-12);
		 cout << "f:   " << k.at(z) + s.at(z) + h.at(z) - f.at(z)(i) << " " << k.at(z) << " " << s.at(z) <<" "<< h.at(z) << endl;
		 cout << "feq: " << k.at(z) + seq.at(z) + heq.at(z) - feq.at(z) << " " << k.at(z) << " " << seq.at(z) <<" "<< heq.at(z) << endl;
		 feq.at(z) = k.at(z) + (2 * seq.at(z) - s.at(z))
		 + ((1 - gamma) * h.at(z) + gamma * heq.at(z));
		 }

		f.at(0)(i) = (1 - beta) * f.at(0)(i) + beta * feq.at(0);
		f.at(1)(i) = (1 - beta) * f.at(1)(i) + beta * feq.at(1);
		f.at(2)(i) = (1 - beta) * f.at(2)(i) + beta * feq.at(2);
		f.at(3)(i) = (1 - beta) * f.at(3)(i) + beta * feq.at(3);
		f.at(4)(i) = (1 - beta) * f.at(4)(i) + beta * feq.at(4);
		f.at(5)(i) = (1 - beta) * f.at(5)(i) + beta * feq.at(5);
		f.at(6)(i) = (1 - beta) * f.at(6)(i) + beta * feq.at(6);
		f.at(7)(i) = (1 - beta) * f.at(7)(i) + beta * feq.at(7);
		f.at(8)(i) = (1 - beta) * f.at(8)(i) + beta * feq.at(8);

	}

}
} /* namespace natrium */
