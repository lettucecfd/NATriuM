/*
 * KBCStandard.cpp
 *
 *  Created on: 17.11.2015
 *      Author: Dominik Wilde
 */

#include "KBCStandard.h"
#include <math.h>
#define KBC_C // if not defined, KBC_D will be used (according to Karlin et al. 2015)
#define EVALUATE_GAMMA // if defined, an  evaluation over time of the stabilizer gamma will be carried out (D2Q9 only)

namespace natrium {

KBCStandard::KBCStandard(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		counter(0), MRT(relaxationParameter, dt, stencil), parameterFile(
				"deviation_KBC_STANDARD.txt") {

}

KBCStandard::~KBCStandard() {

}

double KBCStandard::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

	assert(i < getStencil()->getQ());
	assert(rho > 0);
	assert(u.size() == getStencil()->getD());
	assert(u(0) < 1000000000000000.);
	assert(u(1) < 1000000000000000.);

	double value = getStencil()->getWeight(i) * rho;
	int D = getStencil()->getD();
	numeric_vector c_i = getStencil()->getDirection(i);
	for (int p = 0; p < D; p++) {
		double v = u(p) / (getStencil()->getScaling()) * sqrt(3);

		value *= 2 - sqrt(1 + v * v);
		value *= pow(((2 * v / sqrt(3) + sqrt(1 + v * v)) / (1 - v / sqrt(3))),
				c_i(p) / (getStencil()->getScaling()));

	}
	/*
	 double prefactor = getStencil()->getWeight(i) * rho;
	 double uSquareTerm = -(u * u) / (2 * getStencil()->getSpeedOfSoundSquare());
	 if (0 == i) {
	 cout << value - prefactor * (1 + uSquareTerm);
	 return value;
	 }
	 double mixedTerm = (u * getStencil()->getDirection(i))
	 / getStencil()->getSpeedOfSoundSquare();
	 cout
	 << value
	 - prefactor
	 * (1 + mixedTerm * (1 + 0.5 * (mixedTerm))
	 + uSquareTerm) << endl;
	 */
	return value;

}

void KBCStandard::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
		/*else if (Stencil_D3Q19 == getStencil()->getStencilType()) {
		 collideAllD3Q19(f, densities, velocities, locally_owned_dofs,
		 inInitializationProcedure);*/
	} else if (Stencil_D3Q15 == getStencil()->getStencilType()) {
		collideAllD3Q15(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else {
		throw CollisionException(
				"KBC_Standard only implemented for D2Q9 and D3Q15");
		// Inefficient collision
		//BGK::collideAll(f, densities, velocities, locally_owned_dofs,
		//		inInitializationProcedure);
	}
}

void KBCStandard::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

#define KBC_C // if not defined, KBC_D will be used (according to Karlin et al. 2015)

	size_t Q = getQ();

	stabilizer gamma(locally_owned_dofs.size());
	stabilizer entropy(locally_owned_dofs.size());

	//vector<double> gamma(locally_owned_dofs.size());

	double scaling = getStencil()->getScaling();
//	double cs2 = getStencil()->getSpeedOfSoundSquare() / (scaling * scaling);

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

		if (rho < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

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
//		double scalar_product = ux * ux + uy * uy;

		// moment representation of the populations
		double T = 0, N = 0, Pi_xy = 0, Q_xyy = 0, Q_yxx = 0, A = 0;

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

#ifdef KBC_C

		s.at(0) = -rho * T;
		s.at(1) = 0.5 * rho * 0.5 * (T + N);
		s.at(2) = 0.5 * rho * 0.5 * (T - N);
		s.at(3) = 0.5 * rho * 0.5 * (T + N);
		s.at(4) = 0.5 * rho * 0.5 * (T - N);
		s.at(5) = 0.25 * rho * Pi_xy;
		s.at(6) = 0.25 * rho * -Pi_xy;
		s.at(7) = 0.25 * rho * Pi_xy;
		s.at(8) = 0.25 * rho * -Pi_xy;

#else
		s.at(0) = 0;
		s.at(1) = 0.5 * rho * 0.5 * (+ N);
		s.at(2) = 0.5 * rho * 0.5 * (- N);
		s.at(3) = 0.5 * rho * 0.5 * (+ N);
		s.at(4) = 0.5 * rho * 0.5 * (- N);
		s.at(5) = 0.25 * rho * Pi_xy;
		s.at(6) = 0.25 * rho * -Pi_xy;
		s.at(7) = 0.25 * rho * Pi_xy;
		s.at(8) = 0.25 * rho * -Pi_xy;

#endif

		// higher order part vector h
		vector<double> h(Q);
		/*
		 #ifdef KBC_C
		 h.at(0) = rho * A;
		 h.at(1) = -0.5 * rho * (1 * Q_xyy + A);
		 h.at(3) = -0.5 * rho * (-1 * Q_xyy + A);
		 h.at(2) = -0.5 * rho * (1 * Q_yxx + A);
		 h.at(4) = -0.5 * rho * (-1 * Q_yxx + A);
		 h.at(5) = 0.25 * rho * (1 * Q_xyy + 1 * Q_yxx + A);
		 h.at(6) = 0.25 * rho * (-1 * Q_xyy + 1 * Q_yxx + A);
		 h.at(7) = 0.25 * rho * (-1 * Q_xyy - 1 * Q_yxx + A);
		 h.at(8) = 0.25 * rho * (1 * Q_xyy + -1 * Q_yxx + A);
		 #else
		 h.at(0) = rho * (A - T);
		 h.at(1) = -0.5 * rho * (1 * Q_xyy + A - T);
		 h.at(3) = -0.5 * rho * (-1 * Q_xyy + A - T);
		 h.at(2) = -0.5 * rho * (1 * Q_yxx + A - T);
		 h.at(4) = -0.5 * rho * (-1 * Q_yxx + A - T);
		 h.at(5) = 0.25 * rho * (1 * Q_xyy + 1 * Q_yxx + A);
		 h.at(6) = 0.25 * rho * (-1 * Q_xyy + 1 * Q_yxx + A);
		 h.at(7) = 0.25 * rho * (-1 * Q_xyy - 1 * Q_yxx + A);
		 h.at(8) = 0.25 * rho * (1 * Q_xyy + -1 * Q_yxx + A);
		 #endif
		 */
		h.at(0) = f.at(0)(i) - k.at(0) - s.at(0);
		h.at(1) = f.at(1)(i) - k.at(1) - s.at(1);
		h.at(3) = f.at(3)(i) - k.at(3) - s.at(3);
		h.at(2) = f.at(2)(i) - k.at(2) - s.at(2);
		h.at(4) = f.at(4)(i) - k.at(4) - s.at(4);
		h.at(5) = f.at(5)(i) - k.at(5) - s.at(5);
		h.at(6) = f.at(6)(i) - k.at(6) - s.at(6);
		h.at(7) = f.at(7)(i) - k.at(7) - s.at(7);
		h.at(8) = f.at(8)(i) - k.at(8) - s.at(8);

		// equilibrium vectors for shear vector
		vector<double> seq(Q);
		// equilibrium vectors for higher order part
		vector<double> heq(Q);

		// equilibrium distribution of the population f
		vector<double> feq(Q, 0.0);

		double u_x_i = ux * sqrt(3);
		double u_y_i = uy * sqrt(3);

		double sqrt_ux = sqrt(1 + u_x_i * u_x_i);
		double sqrt_uy = sqrt(1 + u_y_i * u_y_i);

		double prefactor = rho * (2 - sqrt_ux) * (2 - sqrt_uy);

		double postfactor_x = (2 * u_x_i / sqrt(3) + sqrt_ux)
				/ (1 - u_x_i / sqrt(3));
		double postfactor_y = (2 * u_y_i / sqrt(3) + sqrt_uy)
				/ (1 - u_y_i / sqrt(3));

		feq.at(0) = 4. / 9. * prefactor;
		feq.at(1) = 1. / 9. * prefactor * postfactor_x;
		feq.at(2) = 1. / 9. * prefactor * postfactor_y;
		feq.at(3) = 1. / 9. * prefactor / postfactor_x;
		feq.at(4) = 1. / 9. * prefactor / postfactor_y;
		feq.at(5) = 1. / 36. * prefactor * postfactor_x * postfactor_y;
		feq.at(6) = 1. / 36. * prefactor / postfactor_x * postfactor_y;
		feq.at(7) = 1. / 36. * prefactor / postfactor_x / postfactor_y;
		feq.at(8) = 1. / 36. * prefactor * postfactor_x / postfactor_y;

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
#ifdef KBC_C

		seq.at(0) = -rho * T;
		seq.at(1) = 0.5 * rho * 0.5 * (T + N);
		seq.at(3) = 0.5 * rho * 0.5 * (T + N);
		seq.at(2) = 0.5 * rho * 0.5 * (T - N);
		seq.at(4) = 0.5 * rho * 0.5 * (T - N);
		seq.at(5) = 0.25 * rho * Pi_xy;
		seq.at(6) = 0.25 * rho * -Pi_xy;
		seq.at(7) = 0.25 * rho * Pi_xy;
		seq.at(8) = 0.25 * rho * -Pi_xy;

#else
		seq.at(0) = 0;
		seq.at(1) = 0.5 * rho * 0.5 * (+ N);
		seq.at(3) = 0.5 * rho * 0.5 * (+ N);
		seq.at(2) = 0.5 * rho * 0.5 * (- N);
		seq.at(4) = 0.5 * rho * 0.5 * (- N);
		seq.at(5) = 0.25 * rho * Pi_xy;
		seq.at(6) = 0.25 * rho * -Pi_xy;
		seq.at(7) = 0.25 * rho * Pi_xy;
		seq.at(8) = 0.25 * rho * -Pi_xy;
#endif

		/*
		 // calculate higher order equilibrium
		 #ifdef KBC_C

		 heq.at(0) = rho * A;
		 heq.at(1) = -0.5 * rho * (1 * Q_xyy + A);
		 heq.at(3) = -0.5 * rho * (-1 * Q_xyy + A);
		 heq.at(2) = -0.5 * rho * (1 * Q_yxx + A);
		 heq.at(4) = -0.5 * rho * (-1 * Q_yxx + A);
		 heq.at(5) = 0.25 * rho * (1 * Q_xyy + 1 * Q_yxx + A);
		 heq.at(6) = 0.25 * rho * (-1 * Q_xyy + 1 * Q_yxx + A);
		 heq.at(7) = 0.25 * rho * (-1 * Q_xyy - 1 * Q_yxx + A);
		 heq.at(8) = 0.25 * rho * (1 * Q_xyy + -1 * Q_yxx + A);

		 #else
		 heq.at(0) = rho * (A - T);
		 heq.at(1) = -0.5 * rho * (1 * Q_xyy + A - T);
		 heq.at(3) = -0.5 * rho * (-1 * Q_xyy + A - T);
		 heq.at(2) = -0.5 * rho * (1 * Q_yxx + A - T);
		 heq.at(4) = -0.5 * rho * (-1 * Q_yxx + A - T);
		 heq.at(5) = 0.25 * rho * (1 * Q_xyy + 1 * Q_yxx + A);
		 heq.at(6) = 0.25 * rho * (-1 * Q_xyy + 1 * Q_yxx + A);
		 heq.at(7) = 0.25 * rho * (-1 * Q_xyy - 1 * Q_yxx + A);
		 heq.at(8) = 0.25 * rho * (1 * Q_xyy + -1 * Q_yxx + A);
		 #endif
		 */

		heq.at(0) = feq.at(0) - k.at(0) - seq.at(0);
		heq.at(1) = feq.at(1) - k.at(1) - seq.at(1);
		heq.at(3) = feq.at(3) - k.at(3) - seq.at(3);
		heq.at(2) = feq.at(2) - k.at(2) - seq.at(2);
		heq.at(4) = feq.at(4) - k.at(4) - seq.at(4);
		heq.at(5) = feq.at(5) - k.at(5) - seq.at(5);
		heq.at(6) = feq.at(6) - k.at(6) - seq.at(6);
		heq.at(7) = feq.at(7) - k.at(7) - seq.at(7);
		heq.at(8) = feq.at(8) - k.at(8) - seq.at(8);

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

		vector<double> s_pc(Q);
		for (int p = 0; p < Q; p++) {
			s_p	c[p] = s[p]
					- 1. / (getRelaxationParameter() + 0.5) * (s[p] - seq[p]);

			//cout << s_pc[p] << " h: " << h[p] << " " << heq[p] << endl;
		}

		vector<double> psi(Q);
		double sum_0 = 0.0;
		double sum_a = 0.0;
		double sum_b = 0.0;
		double sum_c = 0.0;
		for (int p = 0; p < Q; p++) {
			psi[p] = k[p] + s_pc[p] + (1 - beta) * h[p]
					+ beta * (2 * heq[p] - h[p]);

			sum_0 += (h[p] - heq[p]) * (f.at(p)(i) - psi[p]) / psi[p];

			sum_a += (h[p] - heq[p])
					* (log(psi[p]) - log(getStencil()->getWeight(p))) / psi[p];
			sum_b += (h[p] - heq[p]) * (h[p] - heq[p]) / psi[p];
			sum_c += (h[p] - heq[p]) * (s[p] - seq[p]) / psi[p];
		}
		pout << " sum 0 = " << sum_0 << endl;
		//cout << " sumabc" << sum_a << " " << sum_b << " " << sum_c << endl;

		// stabilizer of KBC model
		//gamma.value.at(i) = 1. / beta - (2 - 1. / beta) * (sum_s / sum_h);
		gamma.value.at(i) = 1. / beta * ((sum_0 + sum_a) / sum_b)
				- 2 * (sum_c / sum_b);
		pout << gamma.value.at(i) << " " << 1. / beta - (2 - 1. / beta) * (sum_s / sum_h) << " ";

		// if the sum_h expression is too small, BGK shall be performed (gamma = 2)
		if (sum_0 < 1e-16 || sum_b < 1e-16 || sum_a < 1e-16 || sum_c < 1e-16) {
			gamma.value.at(i) = 2;

		}

		pout << gamma.value.at(i) << endl;

		entropy.value.at(i) = -(f.at(0)(i) * log(f.at(0)(i) / (4. / 9.))
				+ f.at(1)(i) * log(f.at(1)(i) / (1. / 9.))
				+ f.at(2)(i) * log(f.at(2)(i) / (1. / 9.))
				+ f.at(3)(i) * log(f.at(3)(i) / (1. / 9.))
				+ f.at(4)(i) * log(f.at(4)(i) / (1. / 9.))
				+ f.at(5)(i) * log(f.at(5)(i) / (1. / 36.))
				+ f.at(6)(i) * log(f.at(6)(i) / (1. / 36.))
				+ f.at(7)(i) * log(f.at(7)(i) / (1. / 36.))
				+ f.at(8)(i) * log(f.at(8)(i) / (1. / 36.)));

		// calculate new f
		f.at(0)(i) =
				f.at(0)(i)
						- beta
								* (2 * delta_s.at(0)
										+ gamma.value.at(i) * delta_h.at(0));
		f.at(1)(i) =
				f.at(1)(i)
						- beta
								* (2 * delta_s.at(1)
										+ gamma.value.at(i) * delta_h.at(1));
		f.at(2)(i) =
				f.at(2)(i)
						- beta
								* (2 * delta_s.at(2)
										+ gamma.value.at(i) * delta_h.at(2));
		f.at(3)(i) =
				f.at(3)(i)
						- beta
								* (2 * delta_s.at(3)
										+ gamma.value.at(i) * delta_h.at(3));
		f.at(4)(i) =
				f.at(4)(i)
						- beta
								* (2 * delta_s.at(4)
										+ gamma.value.at(i) * delta_h.at(4));
		f.at(5)(i) =
				f.at(5)(i)
						- beta
								* (2 * delta_s.at(5)
										+ gamma.value.at(i) * delta_h.at(5));
		f.at(6)(i) =
				f.at(6)(i)
						- beta
								* (2 * delta_s.at(6)
										+ gamma.value.at(i) * delta_h.at(6));
		f.at(7)(i) =
				f.at(7)(i)
						- beta
								* (2 * delta_s.at(7)
										+ gamma.value.at(i) * delta_h.at(7));
		f.at(8)(i) =
				f.at(8)(i)
						- beta
								* (2 * delta_s.at(8)
										+ gamma.value.at(i) * delta_h.at(8));

	}

	writeDeviation(gamma.getAverage(), gamma.getDeviation(),
			entropy.getAverage(), entropy.getDeviation());

}

void KBCStandard::collideAllD3Q15(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	// Efficient collision for D2Q9
	size_t Q = 15;
	size_t D = 3;
	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare();
	double prefactor = scaling / cs2;
	double relax_factor = getPrefactor();

	assert(f.size() == Q);
	assert(velocities.size() == D);

#ifdef DEBUG
	size_t n_dofs = f.at(0).size();
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < D; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

	// allocation

	double scalar_product;
	double uSquareTerm;
	double mixedTerm;
	double weighting;
	double rho;
	double ux;
	double uy;
	double uz;
	double T;
	double N_xz;
	double N_yz;
	double Pi_xy;
	double Pi_xz;
	double Pi_yz;
	double Q_xxy, Q_xxz, Q_xyy, Q_xyz, Q_xzz;
	double A;
	double sum_s = 0;
	double sum_h = 0;
	vector<double> feq(Q);
	vector<double> k(Q);
	vector<double> s(Q);
	vector<double> h(Q);
	vector<double> heq(Q);
	vector<double> seq(Q);
	vector<double> delta_s(Q);
	vector<double> delta_h(Q);
	vector<double> delta_seq(Q);
	vector<double> delta_heq(Q);

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// calculate density
		rho = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i)
				+ f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i) + f.at(9)(i)
				+ f.at(10)(i) + f.at(11)(i) + f.at(12)(i) + f.at(13)(i)
				+ f.at(14)(i);
		densities(i) = rho;
		if (rho < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate velocity
			// for all velocity components
			velocities.at(0)(i) = scaling / rho
					* (f.at(1)(i) - f.at(2)(i) + f.at(7)(i) - f.at(8)(i)
							+ f.at(9)(i) - f.at(10)(i) + f.at(11)(i)
							- f.at(12)(i) + f.at(13)(i) - f.at(14)(i));
			velocities.at(1)(i) = scaling / rho
					* (f.at(3)(i) - f.at(4)(i) + f.at(7)(i) - f.at(8)(i)
							+ f.at(9)(i) - f.at(10)(i) - f.at(11)(i)
							+ f.at(12)(i) - f.at(13)(i) + f.at(14)(i));
			velocities.at(2)(i) = scaling / rho
					* (f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i)
							- f.at(9)(i) + f.at(10)(i) + f.at(11)(i)
							- f.at(12)(i) - f.at(13)(i) + f.at(14)(i));
		}
		ux = velocities.at(0)(i) / scaling;
		uy = velocities.at(1)(i) / scaling;
		uz = velocities.at(2)(i) / scaling;

		T = f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i) + f.at(5)(i)
				+ f.at(6)(i) + 3 * f.at(7)(i) + 3 * f.at(8)(i) + 3 * f.at(9)(i)
				+ 3 * f.at(10)(i) + 3 * f.at(11)(i) + 3 * f.at(12)(i)
				+ 3 * f.at(13)(i) + 3 * f.at(14)(i);
		T /= rho;

		//calculate the normal stress difference at unit density (N)
		N_xz = f.at(1)(i) + f.at(2)(i) - f.at(5)(i) - f.at(6)(i);
		N_xz /= rho;

		N_yz = f.at(3)(i) + f.at(4)(i) - f.at(5)(i) - f.at(6)(i);
		N_yz /= rho;

		//calculate the off-diagonal component of the pressure tensor at unit density
		Pi_xy = f.at(7)(i) + f.at(8)(i) + f.at(9)(i) + f.at(10)(i) - f.at(11)(i)
				- f.at(12)(i) - f.at(13)(i) - f.at(14)(i);
		Pi_xy /= rho;

		Pi_xz = f.at(7)(i) + f.at(8)(i) - f.at(9)(i) - f.at(10)(i) + f.at(11)(i)
				+ f.at(12)(i) - f.at(13)(i) - f.at(14)(i);
		Pi_xz /= rho;

		Pi_yz = f.at(7)(i) + f.at(8)(i) - f.at(9)(i) - f.at(10)(i) - f.at(11)(i)
				- f.at(12)(i) + f.at(13)(i) + f.at(14)(i);
		Pi_yz /= rho;

		//calculate moments of third order
		Q_xxy = f.at(7)(i) - f.at(8)(i) + f.at(9)(i) - f.at(10)(i) - f.at(11)(i)
				+ f.at(12)(i) - f.at(13)(i) + f.at(14)(i);
		Q_xxy /= rho;

		Q_xxz = f.at(7)(i) - f.at(8)(i) - f.at(9)(i) + f.at(10)(i) + f.at(11)(i)
				- f.at(12)(i) - f.at(13)(i) + f.at(14)(i);
		Q_xxz /= rho;

		Q_xzz = f.at(7)(i) - f.at(8)(i) + f.at(9)(i) - f.at(10)(i) + f.at(11)(i)
				- f.at(12)(i) + f.at(13)(i) - f.at(14)(i);
		Q_xzz /= rho;

		Q_xyy = f.at(7)(i) - f.at(8)(i) + f.at(9)(i) - f.at(10)(i) + f.at(11)(i)
				- f.at(12)(i) + f.at(13)(i) - f.at(14)(i);
		Q_xyy /= rho;

		Q_xyz = f.at(7)(i) - f.at(8)(i) - f.at(9)(i) + f.at(10)(i) - f.at(11)(i)
				+ f.at(12)(i) + f.at(13)(i) - f.at(14)(i);
		Q_xyz /= rho;

		// calculate moments of forth order
		A = f.at(7)(i) + f.at(8)(i) + f.at(9)(i) + f.at(10)(i) + f.at(11)(i)
				+ f.at(12)(i) + f.at(13)(i) + f.at(14)(i);
		A = A / rho;

		// kinematic part vector k
		k.at(0) = rho;
		k.at(1) = rho / 6 * (3 * ux);
		k.at(2) = rho / 6 * (3 * -ux);
		k.at(3) = rho / 6 * (3 * uy);
		k.at(4) = rho / 6 * (3 * -uy);
		k.at(5) = rho / 6 * (3 * uz);
		k.at(6) = rho / 6 * (3 * -uz);
		k.at(7) = 0;
		k.at(8) = 0;
		k.at(9) = 0;
		k.at(10) = 0;
		k.at(11) = 0;
		k.at(12) = 0;
		k.at(13) = 0;
		k.at(14) = 0;

		// shear part vector s
		s.at(0) = rho * -T;
		s.at(1) = 1. / 6. * rho * (2 * N_xz - N_yz + T);
		s.at(2) = s.at(1);
		s.at(3) = 1. / 6. * rho * (-N_xz + 2 * N_yz + T);
		s.at(4) = s.at(3);
		s.at(5) = 1. / 6. * rho * (-N_xz - N_yz + T);
		s.at(6) = s.at(5);
		s.at(7) = 1. / 8. * rho * Q_xyz;
		s.at(8) = -s.at(7);
		s.at(9) = -s.at(7);
		s.at(10) = s.at(7);
		s.at(11) = -s.at(7);
		s.at(12) = s.at(7);
		s.at(13) = s.at(7);
		s.at(14) = -s.at(7);

		/*		// higher order part vector h
		 h.at(0) = 2. * A * rho;
		 h.at(1) = 1. / 6. * rho * 3 * (-Q_xyy - A);
		 h.at(2) = 1. / 6. * rho * 3 * (Q_xyy - A);
		 h.at(3) = 1. / 6. * rho * 3 * (-Q_xxy - A);
		 h.at(4) = 1. / 6. * rho * 3 * (Q_xxy - A);
		 h.at(5) = 1. / 6. * rho * 3 * (-Q_xxz - A);
		 h.at(6) = 1. / 6. * rho * 3 * (Q_xxz - A);
		 h.at(7) = 1. / 8. * rho
		 * (Q_xzz + Q_xxy + Q_xxz + Pi_xy + Pi_xz + Pi_yz + A);
		 h.at(8) = 1. / 8. * rho
		 * (-Q_xzz - Q_xxy - Q_xxz + Pi_xy + Pi_xz + Pi_yz + A);

		 h.at(9) = 1. / 8. * rho
		 * (1 * Q_xzz + 1 * Q_xxy + -1 * Q_xxz + 1 * 1 * Pi_xy
		 + 1 * -1 * Pi_xz + 1 * -1 * Pi_yz + A);
		 h.at(10) = 1. / 8. * rho
		 * (-1 * Q_xzz + -1 * Q_xxy + 1 * Q_xxz + -1 * -1 * Pi_xy
		 + -1 * 1 * Pi_xz + -1 * 1 * Pi_yz + A);
		 h.at(11) = 1. / 8. * rho
		 * (1 * Q_xzz + -1 * Q_xxy + 1 * Q_xxz + 1 * -1 * Pi_xy
		 + 1 * 1 * Pi_xz + -1 * 1 * Pi_yz + A);
		 h.at(12) = 1. / 8. * rho
		 * (-1 * Q_xzz + 1 * Q_xxy + -1 * Q_xxz + -1 * 1 * Pi_xy
		 + -1 * -1 * Pi_xz + 1 * -1 * Pi_yz + A);
		 h.at(13) = 1. / 8. * rho
		 * (1 * Q_xzz + -1 * Q_xxy + -1 * Q_xxz + 1 * -1 * Pi_xy
		 + 1 * -1 * Pi_xz + -1 * -1 * Pi_yz + A);

		 h.at(14) = 1. / 8. * rho
		 * (-1 * Q_xzz + 1 * Q_xxy + 1 * Q_xxz + -1 * 1 * Pi_xy
		 + -1 * 1 * Pi_xz + 1 * 1 * Pi_yz + A);*/

		h.at(0) = f.at(0)(i) - k.at(0) - s.at(0);
		h.at(1) = f.at(1)(i) - k.at(1) - s.at(1);
		h.at(2) = f.at(2)(i) - k.at(2) - s.at(2);
		h.at(3) = f.at(3)(i) - k.at(3) - s.at(3);
		h.at(4) = f.at(4)(i) - k.at(4) - s.at(4);
		h.at(5) = f.at(5)(i) - k.at(5) - s.at(5);
		h.at(6) = f.at(6)(i) - k.at(6) - s.at(6);
		h.at(7) = f.at(7)(i) - k.at(7) - s.at(7);
		h.at(8) = f.at(8)(i) - k.at(8) - s.at(8);
		h.at(9) = f.at(9)(i) - k.at(9) - s.at(9);
		h.at(10) = f.at(10)(i) - k.at(10) - s.at(10);
		h.at(11) = f.at(11)(i) - k.at(11) - s.at(11);
		h.at(12) = f.at(12)(i) - k.at(12) - s.at(12);
		h.at(13) = f.at(13)(i) - k.at(13) - s.at(13);
		h.at(14) = f.at(14)(i) - k.at(14) - s.at(14);

		// calculate equilibrium distribution
		scalar_product = ux * ux + uy * uy + uz * uz;
		uSquareTerm = -scalar_product / (2 * cs2);
		// direction 0
		weighting = 2. / 9. * rho;
		feq.at(0) = weighting * (1 + uSquareTerm);
		// directions 1-6 (The mixed term is (e_i *u)^2 / (2c_s^4)
		weighting = 1. / 9. * rho;
		mixedTerm = prefactor * (velocities.at(0)(i));
		feq.at(1) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(2) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * ((velocities.at(1)(i)));
		feq.at(3) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(4) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * ((velocities.at(2)(i)));
		feq.at(5) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(6) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		// directions 7-18
		weighting = 1. / 72. * rho;
		mixedTerm = prefactor
				* ((velocities.at(0)(i)) + (velocities.at(1)(i))
						+ (velocities.at(2)(i)));
		feq.at(7) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(8) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor
				* ((velocities.at(0)(i)) + (velocities.at(1)(i))
						- (velocities.at(2)(i)));
		feq.at(9) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(10) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor
				* ((velocities.at(0)(i)) - (velocities.at(1)(i))
						+ (velocities.at(2)(i)));
		feq.at(11) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(12) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor
				* ((velocities.at(0)(i)) - (velocities.at(1)(i))
						- (velocities.at(2)(i)));
		feq.at(13) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(14) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);

		// calculate all moments for the equilibrium distribution function again

		//calculate the trace of the pressure tensor at unit density (T)
		T = feq.at(1) + feq.at(2) + feq.at(3) + feq.at(4) + feq.at(5)
				+ feq.at(6) + 3 * feq.at(7) + 3 * feq.at(8) + 3 * feq.at(9)
				+ 3 * feq.at(10) + 3 * feq.at(11) + 3 * feq.at(12)
				+ 3 * feq.at(13) + 3 * feq.at(14);
		T /= rho;

		//calculate the normal stress difference at unit density (N)
		N_xz = feq.at(1) + feq.at(2) - feq.at(5) - feq.at(6);
		N_xz /= rho;

		N_yz = feq.at(3) + feq.at(4) - feq.at(5) - feq.at(6);
		N_yz /= rho;

		//calculate the off-diagonal component of the pressure tensor at unit density
		Pi_xy = feq.at(7) + feq.at(8) + feq.at(9) + feq.at(10) - feq.at(11)
				- feq.at(12) - feq.at(13) - feq.at(14);
		Pi_xy /= rho;

		Pi_xz = feq.at(7) + feq.at(8) - feq.at(9) - feq.at(10) + feq.at(11)
				+ feq.at(12) - feq.at(13) - feq.at(14);
		Pi_xz /= rho;

		Pi_yz = feq.at(7) + feq.at(8) - feq.at(9) - feq.at(10) - feq.at(11)
				- feq.at(12) + feq.at(13) + feq.at(14);
		Pi_yz /= rho;

		//calculate moments of third order
		Q_xxy = feq.at(7) - feq.at(8) + feq.at(9) - feq.at(10) - feq.at(11)
				+ feq.at(12) - feq.at(13) + feq.at(14);
		Q_xxy /= rho;

		Q_xxz = feq.at(7) - feq.at(8) - feq.at(9) + feq.at(10) + feq.at(11)
				- feq.at(12) - feq.at(13) + feq.at(14);
		Q_xxz /= rho;

		Q_xzz = feq.at(7) - feq.at(8) + feq.at(9) - feq.at(10) + feq.at(11)
				- feq.at(12) + feq.at(13) - feq.at(14);
		Q_xzz /= rho;

		Q_xyy = feq.at(7) - feq.at(8) + feq.at(9) - feq.at(10) + feq.at(11)
				- feq.at(12) + feq.at(13) - feq.at(14);
		Q_xyy /= rho;

		Q_xyz = feq.at(7) - feq.at(8) - feq.at(9) + feq.at(10) - feq.at(11)
				+ feq.at(12) + feq.at(13) - feq.at(14);
		Q_xyz /= rho;

		A = feq.at(7) + feq.at(8) + feq.at(9) + feq.at(10) + feq.at(11)
				+ feq.at(12) + feq.at(13) + feq.at(14);
		A = A / rho;

		// calculate shear equilibrium
		seq.at(0) = rho * -T;
		seq.at(1) = 1. / 6. * rho * (2 * N_xz - N_yz + T);
		seq.at(2) = seq.at(1);
		seq.at(3) = 1. / 6. * rho * (-N_xz + 2 * N_yz + T);
		seq.at(4) = seq.at(3);
		seq.at(5) = 1. / 6. * rho * (-N_xz - N_yz + T);
		seq.at(6) = seq.at(5);
		seq.at(7) = 1. / 8. * rho * Q_xyz;
		seq.at(8) = -seq.at(7);
		seq.at(9) = -seq.at(7);
		seq.at(10) = seq.at(7);
		seq.at(11) = -seq.at(7);
		seq.at(12) = seq.at(7);
		seq.at(13) = seq.at(7);
		seq.at(14) = -seq.at(7);

		/*		// calculate higher order equilibrium
		 heq.at(0) = 2. * A * rho;
		 heq.at(1) = 1. / 6. * rho * 3 * (-Q_xyy - A);
		 heq.at(2) = 1. / 6. * rho * 3 * (Q_xyy - A);
		 heq.at(3) = 1. / 6. * rho * 3 * (-Q_xxy - A);
		 heq.at(4) = 1. / 6. * rho * 3 * (Q_xxy - A);
		 heq.at(5) = 1. / 6. * rho * 3 * (-Q_xxz - A);
		 heq.at(6) = 1. / 6. * rho * 3 * (Q_xxz - A);
		 heq.at(7) = 1. / 8. * rho
		 * (Q_xzz + Q_xxy + Q_xxz + Pi_xy + Pi_xz + Pi_yz + A);
		 heq.at(8) = 1. / 8. * rho
		 * (-Q_xzz - Q_xxy - Q_xxz + Pi_xy + Pi_xz + Pi_yz + A);

		 heq.at(9) = 1. / 8. * rho
		 * (1 * Q_xzz + 1 * Q_xxy + -1 * Q_xxz + 1 * 1 * Pi_xy
		 + 1 * -1 * Pi_xz + 1 * -1 * Pi_yz + A);
		 heq.at(10) = 1. / 8. * rho
		 * (-1 * Q_xzz + -1 * Q_xxy + 1 * Q_xxz + -1 * -1 * Pi_xy
		 + -1 * 1 * Pi_xz + -1 * 1 * Pi_yz + A);
		 heq.at(11) = 1. / 8. * rho
		 * (1 * Q_xzz + -1 * Q_xxy + 1 * Q_xxz + 1 * -1 * Pi_xy
		 + 1 * 1 * Pi_xz + -1 * 1 * Pi_yz + A);
		 heq.at(12) = 1. / 8. * rho
		 * (-1 * Q_xzz + 1 * Q_xxy + -1 * Q_xxz + -1 * 1 * Pi_xy
		 + -1 * -1 * Pi_xz + 1 * -1 * Pi_yz + A);
		 heq.at(13) = 1. / 8. * rho
		 * (1 * Q_xzz + -1 * Q_xxy + -1 * Q_xxz + 1 * -1 * Pi_xy
		 + 1 * -1 * Pi_xz + -1 * -1 * Pi_yz + A);

		 heq.at(14) = 1. / 8. * rho
		 * (-1 * Q_xzz + 1 * Q_xxy + 1 * Q_xxz + -1 * 1 * Pi_xy
		 + -1 * 1 * Pi_xz + 1 * 1 * Pi_yz + A);*/

		heq.at(0) = feq.at(0) - k.at(0) - seq.at(0);
		heq.at(1) = feq.at(1) - k.at(1) - seq.at(1);
		heq.at(2) = feq.at(2) - k.at(2) - seq.at(2);
		heq.at(3) = feq.at(3) - k.at(3) - seq.at(3);
		heq.at(4) = feq.at(4) - k.at(4) - seq.at(4);
		heq.at(5) = feq.at(5) - k.at(5) - seq.at(5);
		heq.at(6) = feq.at(6) - k.at(6) - seq.at(6);
		heq.at(7) = feq.at(7) - k.at(7) - seq.at(7);
		heq.at(8) = feq.at(8) - k.at(8) - seq.at(8);
		heq.at(9) = feq.at(9) - k.at(9) - seq.at(9);
		heq.at(10) = feq.at(10) - k.at(10) - seq.at(10);
		heq.at(11) = feq.at(11) - k.at(11) - seq.at(11);
		heq.at(12) = feq.at(12) - k.at(12) - seq.at(12);
		heq.at(13) = feq.at(13) - k.at(13) - seq.at(13);
		heq.at(14) = feq.at(14) - k.at(14) - seq.at(14);

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
		delta_s.at(9) = s.at(9) - seq.at(9);
		delta_s.at(10) = s.at(10) - seq.at(10);
		delta_s.at(11) = s.at(11) - seq.at(11);
		delta_s.at(12) = s.at(12) - seq.at(12);
		delta_s.at(13) = s.at(13) - seq.at(13);
		delta_s.at(14) = s.at(14) - seq.at(14);

		delta_h.at(0) = h.at(0) - heq.at(0);
		delta_h.at(1) = h.at(1) - heq.at(1);
		delta_h.at(2) = h.at(2) - heq.at(2);
		delta_h.at(3) = h.at(3) - heq.at(3);
		delta_h.at(4) = h.at(4) - heq.at(4);
		delta_h.at(5) = h.at(5) - heq.at(5);
		delta_h.at(6) = h.at(6) - heq.at(6);
		delta_h.at(7) = h.at(7) - heq.at(7);
		delta_h.at(8) = h.at(8) - heq.at(8);
		delta_h.at(9) = h.at(9) - heq.at(9);
		delta_h.at(10) = h.at(10) - heq.at(10);
		delta_h.at(11) = h.at(11) - heq.at(11);
		delta_h.at(12) = h.at(12) - heq.at(12);
		delta_h.at(13) = h.at(13) - heq.at(13);
		delta_h.at(14) = h.at(14) - heq.at(14);

		delta_seq.at(0) = delta_s.at(0) * delta_h.at(0) / feq.at(0);
		delta_seq.at(1) = delta_s.at(1) * delta_h.at(1) / feq.at(1);
		delta_seq.at(2) = delta_s.at(2) * delta_h.at(2) / feq.at(2);
		delta_seq.at(3) = delta_s.at(3) * delta_h.at(3) / feq.at(3);
		delta_seq.at(4) = delta_s.at(4) * delta_h.at(4) / feq.at(4);
		delta_seq.at(5) = delta_s.at(5) * delta_h.at(5) / feq.at(5);
		delta_seq.at(6) = delta_s.at(6) * delta_h.at(6) / feq.at(6);
		delta_seq.at(7) = delta_s.at(7) * delta_h.at(7) / feq.at(7);
		delta_seq.at(8) = delta_s.at(8) * delta_h.at(8) / feq.at(8);
		delta_seq.at(9) = delta_s.at(9) * delta_h.at(9) / feq.at(9);
		delta_seq.at(10) = delta_s.at(10) * delta_h.at(10) / feq.at(10);
		delta_seq.at(11) = delta_s.at(11) * delta_h.at(11) / feq.at(11);
		delta_seq.at(12) = delta_s.at(12) * delta_h.at(12) / feq.at(12);
		delta_seq.at(13) = delta_s.at(13) * delta_h.at(13) / feq.at(13);
		delta_seq.at(14) = delta_s.at(14) * delta_h.at(14) / feq.at(14);

		delta_heq.at(0) = delta_h.at(0) * delta_h.at(0) / feq.at(0);
		delta_heq.at(1) = delta_h.at(1) * delta_h.at(1) / feq.at(1);
		delta_heq.at(2) = delta_h.at(2) * delta_h.at(2) / feq.at(2);
		delta_heq.at(3) = delta_h.at(3) * delta_h.at(3) / feq.at(3);
		delta_heq.at(4) = delta_h.at(4) * delta_h.at(4) / feq.at(4);
		delta_heq.at(5) = delta_h.at(5) * delta_h.at(5) / feq.at(5);
		delta_heq.at(6) = delta_h.at(6) * delta_h.at(6) / feq.at(6);
		delta_heq.at(7) = delta_h.at(7) * delta_h.at(7) / feq.at(7);
		delta_heq.at(8) = delta_h.at(8) * delta_h.at(8) / feq.at(8);
		delta_heq.at(9) = delta_h.at(9) * delta_h.at(9) / feq.at(9);
		delta_heq.at(10) = delta_h.at(10) * delta_h.at(10) / feq.at(10);
		delta_heq.at(11) = delta_h.at(11) * delta_h.at(11) / feq.at(11);
		delta_heq.at(12) = delta_h.at(12) * delta_h.at(12) / feq.at(12);
		delta_heq.at(13) = delta_h.at(13) * delta_h.at(13) / feq.at(13);
		delta_heq.at(14) = delta_h.at(14) * delta_h.at(14) / feq.at(14);

		// sum of all entropic scalar products of s
		sum_s = delta_seq.at(0) + delta_seq.at(1) + delta_seq.at(2)
				+ delta_seq.at(3) + delta_seq.at(4) + delta_seq.at(5)
				+ delta_seq.at(6) + delta_seq.at(7) + delta_seq.at(8)
				+ delta_seq.at(9) + delta_seq.at(10) + delta_seq.at(11)
				+ delta_seq.at(12) + delta_seq.at(13) + delta_seq.at(14);

		// sum of all entropic scalar products of h
		sum_h = delta_heq.at(0) + delta_heq.at(1) + delta_heq.at(2)
				+ delta_heq.at(3) + delta_heq.at(4) + delta_heq.at(5)
				+ delta_heq.at(6) + delta_heq.at(7) + delta_heq.at(8)
				+ delta_heq.at(9) + delta_heq.at(10) + delta_heq.at(11)
				+ delta_heq.at(12) + delta_heq.at(13) + delta_heq.at(14);

		// relaxation parameter (2*beta = omega of BGK)
		double beta = 1. / (getRelaxationParameter() + 0.5) / 2;

		// stabilizer of KBC model
		double gamma = 1. / beta - (2 - 1. / beta) * sum_s / sum_h;

		// if the sum_h expression is too small, BGK shall be performed (gamma = 2)
		if (sum_h < 1e-16) {
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
		f.at(9)(i) = f.at(9)(i)
				- beta * (2 * delta_s.at(9) + gamma * delta_h.at(9));
		f.at(10)(i) = f.at(10)(i)
				- beta * (2 * delta_s.at(10) + gamma * delta_h.at(10));
		f.at(11)(i) = f.at(11)(i)
				- beta * (2 * delta_s.at(11) + gamma * delta_h.at(11));
		f.at(12)(i) = f.at(12)(i)
				- beta * (2 * delta_s.at(12) + gamma * delta_h.at(12));
		f.at(13)(i) = f.at(13)(i)
				- beta * (2 * delta_s.at(13) + gamma * delta_h.at(13));
		f.at(14)(i) = f.at(14)(i)
				- beta * (2 * delta_s.at(14) + gamma * delta_h.at(14));

	}
}
}
/* namespace natrium */
