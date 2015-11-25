/*
 * KBCStandard.cpp
 *
 *  Created on: 17.11.2015
 *      Author: dominik
 */

#include "KBCStandard.h"

namespace natrium {

KBCStandard::KBCStandard(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil) {

}

KBCStandard::~KBCStandard() {

}

void KBCStandard::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		bool inInitializationProcedure = false) const {

	if (Stencil_D2Q9 != getStencil()->getStencilType()) {
		// No MRT collision for other StencilTypes than D2Q9
		throw CollisionException("MRT only implemented for D2Q9");
	} else {

		double cs2 = getStencil()->getSpeedOfSoundSquare();

		size_t n_dofs = f.at(0).size();

		size_t Q = getQ();

		double scaling = getStencil()->getScaling();

		for (size_t i = 0; i < n_dofs; i++) {

			// calculate density
			densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);


			if (not inInitializationProcedure) {

				velocities.at(0)(i) = scaling / densities(i)
						* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
								- f.at(6)(i) - f.at(7)(i));
				velocities.at(1)(i) = scaling / densities(i)
						* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
								- f.at(7)(i) - f.at(8)(i));
			}

			double scalar_product = velocities.at(0)(i) * velocities.at(0)(i)
					+ velocities.at(1)(i) * velocities.at(1)(i);

			double T = 0, N = 0, Pi_xy = 0, Q_xyy = 0, Q_xxy = 0, Q_yxx = 0, A =
					0;

			T = f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i)
					+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));

			T = T / densities(i);


			N = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);

			N = N / densities(i);


			Pi_xy = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);
			Pi_xy = Pi_xy / densities(i);

			Q_yxx = f.at(5)(i) + f.at(6)(i) - f.at(7)(i) - f.at(8)(i);
			Q_yxx = Q_yxx / densities(i);

			A = f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
			A = A / densities(i);

			T = T - scalar_product;


			N = N
					- (velocities.at(0)(i) * velocities.at(0)(i)
							- velocities.at(1)(i) * velocities.at(1)(i));


			Pi_xy = Pi_xy - velocities.at(0)(i) * velocities.at(1)(i);



			Q_xyy = Q_xyy - 2 * velocities.at(1)(i) * Pi_xy
					+ 0.5 * velocities.at(0)(i) * N
					- 0.5 * velocities.at(0)(i) * T
					- velocities.at(0)(i) * velocities.at(1)(i)
							* velocities.at(1)(i);
			Q_yxx = Q_yxx - 2 * velocities.at(0)(i) * Pi_xy
					- 0.5 * velocities.at(1)(i) * N
					- 0.5 * velocities.at(1)(i) * T
					- velocities.at(1)(i) * velocities.at(0)(i)
							* velocities.at(0)(i);



			A = A
					- 2
							* (velocities.at(0)(i) * Q_xyy
									+ velocities.at(1)(i) * Q_yxx)
					- 4 * velocities.at(0)(i) * velocities.at(1)(i) * Pi_xy
					- 0.5 * scalar_product * T
					+ 0.5
							* (velocities.at(0)(i) * velocities.at(0)(i)
									- velocities.at(1)(i) * velocities.at(1)(i))
							* N
					- velocities.at(0)(i) * velocities.at(0)(i)
							* velocities.at(1)(i) * velocities.at(1)(i);



			vector<double> k(Q);

			k.at(0) = densities(i) * (1 - scalar_product);
			k.at(1) = 0.5 * densities(i)
					* (velocities.at(0)(i) * velocities.at(0)(i)
							+ 1 * velocities.at(0)(i));
			k.at(3) = 0.5 * densities(i)
					* (velocities.at(0)(i) * velocities.at(0)(i)
							- 1 * velocities.at(0)(i));
			k.at(2) = 0.5 * densities(i)
					* (velocities.at(1)(i) * velocities.at(1)(i)
							+ 1 * velocities.at(1)(i));
			k.at(4) = 0.5 * densities(i)
					* (velocities.at(1)(i) * velocities.at(1)(i)
							- 1 * velocities.at(1)(i));
			k.at(5) = 0.25 * densities(i) * (1 * 1) * velocities.at(0)(i)
					* velocities.at(1)(i);
			k.at(6) = 0.25 * densities(i) * (1 * -1) * velocities.at(0)(i)
					* velocities.at(1)(i);
			k.at(7) = 0.25 * densities(i) * (-1 * -1) * velocities.at(0)(i)
					* velocities.at(1)(i);
			k.at(8) = 0.25 * densities(i) * (-1 * 1) * velocities.at(0)(i)
					* velocities.at(1)(i);





			vector<double> s(Q);
#define RHO densities(i)
			s.at(0) = RHO
					* ((4 * velocities.at(0)(i)) * velocities.at(1)(i) * Pi_xy
							- 0.5
									* (velocities.at(0)(i) * velocities.at(0)(i)
											- velocities.at(1)(i)
													* velocities.at(1)(i)) * N
							+ 0.5 * (scalar_product - 2) * T);
			s.at(1) = 0.5 * RHO
					* (0.5
							* (1 + 1 * velocities.at(0)(i)
									+ velocities.at(0)(i) * velocities.at(0)(i)
									- velocities.at(1)(i) * velocities.at(1)(i))
							* N
							- (2 * 1 * velocities.at(1)(i)
									+ 4 * velocities.at(0)(i)
											* velocities.at(1)(i)) * Pi_xy
							+ 0.5
									* (1 - 1 * velocities.at(0)(i)
											- scalar_product) * T);

			s.at(3) = 0.5 * RHO
					* (0.5
							* (1 - 1 * velocities.at(0)(i)
									+ velocities.at(0)(i) * velocities.at(0)(i)
									- velocities.at(1)(i) * velocities.at(1)(i))
							* N
							- (2 * -1 * velocities.at(1)(i)
									+ 4 * velocities.at(0)(i)
											* velocities.at(1)(i)) * Pi_xy
							+ 0.5
									* (1 + 1 * velocities.at(0)(i)
											- scalar_product) * T);

			s.at(2) = 0.5 * RHO
					* (0.5
							* (-1 - 1 * velocities.at(1)(i)
									+ velocities.at(0)(i) * velocities.at(0)(i)
									- velocities.at(1)(i) * velocities.at(1)(i))
							* N
							- (2 * 1 * velocities.at(0)(i)
									+ 4 * velocities.at(0)(i)
											* velocities.at(1)(i)) * Pi_xy
							+ 0.5
									* (1 - 1 * velocities.at(1)(i)
											- scalar_product) * T);

			s.at(4) = 0.5 * RHO
					* (0.5
							* (-1 + 1 * velocities.at(1)(i)
									+ velocities.at(0)(i) * velocities.at(0)(i)
									- velocities.at(1)(i) * velocities.at(1)(i))
							* N
							- (2 * -1 * velocities.at(0)(i)
									+ 4 * velocities.at(0)(i)
											* velocities.at(1)(i)) * Pi_xy
							+ 0.5
									* (1 + 1 * velocities.at(1)(i)
											- scalar_product) * T);

			s.at(5) = 0.25 * RHO
					* ((4 * velocities.at(0)(i) * velocities.at(1)(i)
							+ (1) * (1) + 2 * 1 * velocities.at(1)(i)
							+ 2 * 1 * velocities.at(0)(i)) * Pi_xy
							+ 0.5
									* (-velocities.at(0)(i)
											* velocities.at(0)(i)
											+ velocities.at(1)(i)
													* velocities.at(1)(i)
											- 1 * velocities.at(0)(i)
											+ 1 * velocities.at(1)(i)) * N
							+ 0.5
									* (scalar_product + 1 * velocities.at(0)(i)
											+ 1 * velocities.at(1)(i)) * T);
			s.at(6) = 0.25 * RHO
					* ((4 * velocities.at(0)(i) * velocities.at(1)(i)
							+ (-1) * (1) + 2 * -1 * velocities.at(1)(i)
							+ 2 * 1 * velocities.at(0)(i)) * Pi_xy
							+ 0.5
									* (-velocities.at(0)(i)
											* velocities.at(0)(i)
											+ velocities.at(1)(i)
													* velocities.at(1)(i)
											+ 1 * velocities.at(0)(i)
											+ 1 * velocities.at(1)(i)) * N
							+ 0.5
									* (scalar_product
											+ (-1) * velocities.at(0)(i)
											+ 1 * velocities.at(1)(i)) * T);
			s.at(7) = 0.25 * RHO
					* ((4 * velocities.at(0)(i) * velocities.at(1)(i)
							+ (-1) * (-1) + 2 * -1 * velocities.at(1)(i)
							+ 2 * (-1) * velocities.at(0)(i)) * Pi_xy
							+ 0.5
									* (-velocities.at(0)(i)
											* velocities.at(0)(i)
											+ velocities.at(1)(i)
													* velocities.at(1)(i)
											+ 1 * velocities.at(0)(i)
											- 1 * velocities.at(1)(i)) * N
							+ 0.5
									* (scalar_product
											+ (-1) * velocities.at(0)(i)
											- 1 * velocities.at(1)(i)) * T);

			s.at(8) = 0.25 * RHO
					* ((4 * velocities.at(0)(i) * velocities.at(1)(i)
							+ (1) * (-1) + 2 * 1 * velocities.at(1)(i)
							+ 2 * (-1) * velocities.at(0)(i)) * Pi_xy
							+ 0.5
									* (-velocities.at(0)(i)
											* velocities.at(0)(i)
											+ velocities.at(1)(i)
													* velocities.at(1)(i)
											- 1 * velocities.at(0)(i)
											- 1 * velocities.at(1)(i)) * N
							+ 0.5
									* (scalar_product
											+ (1) * velocities.at(0)(i)
											- 1 * velocities.at(1)(i)) * T);





			vector<double> h(Q);

			h.at(0) = RHO
					* (2 * velocities.at(0)(i) * Q_xyy
							+ 2 * velocities.at(1)(i) * Q_yxx + A);
			h.at(1) = 0.5 * RHO
					* (-(1 + 2 * velocities.at(0)(i)) * Q_xyy
							- 2 * velocities.at(1)(i) * Q_yxx - A);
			h.at(3) = 0.5 * RHO
					* (-(-1 + 2 * velocities.at(0)(i)) * Q_xyy
							- 2 * velocities.at(1)(i) * Q_yxx - A);
			h.at(2) = 0.5 * RHO
					* (-(1 + 2 * velocities.at(1)(i)) * Q_yxx
							- 2 * velocities.at(0)(i) * Q_xyy - A);
			h.at(4) = 0.5 * RHO
					* (-(-1 + 2 * velocities.at(1)(i)) * Q_yxx
							- 2 * velocities.at(0)(i) * Q_xyy - A);
			h.at(5) = 0.25 * RHO
					* ((1 + 2 * velocities.at(0)(i)) * Q_xyy
							+ (1 + 2 * velocities.at(1)(i)) * Q_yxx + A);
			h.at(6) = 0.25 * RHO
					* ((-1 + 2 * velocities.at(0)(i)) * Q_xyy
							+ (1 + 2 * velocities.at(1)(i)) * Q_yxx + A);
			h.at(7) = 0.25 * RHO
					* ((-1 + 2 * velocities.at(0)(i)) * Q_xyy
							+ (-1 + 2 * velocities.at(1)(i)) * Q_yxx + A);
			h.at(8) = 0.25 * RHO
					* ((1 + 2 * velocities.at(0)(i)) * Q_xyy
							+ (-1 + 2 * velocities.at(1)(i)) * Q_yxx + A);


			vector<double> seq(Q);
			vector<double> heq(Q);

			T = 2 * cs2;
			A = cs2 * cs2;

			seq.at(0) = RHO * (0.5 * (scalar_product - 2) * T);
			seq.at(1) = 0.5 * RHO * 0.5
					* (1 - velocities.at(0)(i) - scalar_product) * T;
			seq.at(3) = 0.5 * RHO * 0.5
					* (1 + velocities.at(0)(i) - scalar_product) * T;
			seq.at(2) = 0.5 * RHO * 0.5
					* (1 - velocities.at(1)(i) - scalar_product) * T;
			seq.at(4) = 0.5 * RHO * 0.5
					* (1 + velocities.at(1)(i) - scalar_product) * T;
			seq.at(5) = 0.25 * RHO * 0.5
					* (scalar_product + velocities.at(0)(i)
							+ velocities.at(1)(i)) * T;
			seq.at(6) = 0.25 * RHO * 0.5
					* (scalar_product - velocities.at(0)(i)
							+ velocities.at(1)(i)) * T;
			seq.at(7) = 0.25 * RHO * 0.5
					* (scalar_product - velocities.at(0)(i)
							- velocities.at(1)(i)) * T;
			seq.at(8) = 0.25 * RHO * 0.5
					* (scalar_product + velocities.at(0)(i)
							- velocities.at(1)(i)) * T;

			heq.at(0) = RHO * A;
			heq.at(1) = heq.at(3) = -0.5 * RHO * A;
			heq.at(2) = heq.at(4) = -0.5 * RHO * A;
			heq.at(5) = heq.at(6) = heq.at(7) = heq.at(8) = 0.25 * RHO * A;



			vector<double> feq(Q, 0.0);
			double uSquareTerm;
			double mixedTerm;
			double weighting;
			double prefactor = scaling / cs2;

			uSquareTerm = -scalar_product / (2 * cs2);
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
			mixedTerm = prefactor
					* (-velocities.at(0)(i) + velocities.at(1)(i));
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
					+ delta_h.at(3) + delta_h.at(4) + delta_h.at(5)
					+ delta_h.at(6) + delta_h.at(7) + delta_h.at(8);

			double sum_s = delta_s.at(0) + delta_s.at(1) + delta_s.at(2)
					+ delta_s.at(3) + delta_s.at(4) + delta_s.at(5)
					+ delta_s.at(6) + delta_s.at(7) + delta_s.at(8);

			double viscosity = getRelaxationParameter()*cs2;

			double beta = 1. / (viscosity / cs2 + 0.5)/2;

			double gamma = 1. / beta - (2 - 1 / beta) * sum_s / sum_h;

			if (sum_h < 0.000001)
			{
				gamma = 1. / beta - (2 - 1 / beta);
			}

			for (int z = 0; z < 9; z++) {
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
}
}

/* namespace natrium */
