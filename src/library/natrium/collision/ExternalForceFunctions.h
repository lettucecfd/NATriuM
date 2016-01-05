/**
 * @file ExternalForceFunctions.h
 * @short Inline functions to incorporate external forces
 * @date 04.01.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef EXTERNALFORCEFUNCTIONS_H_
#define EXTERNALFORCEFUNCTIONS_H_

namespace natrium {

namespace ExternalForceFunctions {

inline void applyExactDifferenceForcingD2Q9(double f_i[9], double force_x, double force_y,
		double u_0_i, double u_1_i, double rho_i, double dt, double cs2,
		double prefactor) {
	double mixed_term_shifted;
	double mixed_term_eq;
	double s_i[9];

	// Calculate velocity shift
	double deltaU_0 = force_x * dt / rho_i;
	double deltaU_1 = force_y * dt / rho_i;

	double u_eq_shifted0 = u_0_i + deltaU_0;
	double u_eq_shifted1 = u_1_i + deltaU_1;

	double scalar_product_eq_shifted = u_eq_shifted0 * u_eq_shifted0
			+ u_eq_shifted1 * u_eq_shifted1;
	double uSquareTerm_eq_shifted = -scalar_product_eq_shifted / (2 * cs2);

	double scalar_product_eq = u_0_i * u_0_i + u_1_i * u_1_i;
	double u_square_term_eq = -scalar_product_eq / (2 * cs2);

	// Calculate source term
	// direction 0
	double weighting = 4. / 9. * rho_i;
	s_i[0] = weighting * (1 + uSquareTerm_eq_shifted)
			- weighting * (1 + u_square_term_eq);

	// directions 1-4
	weighting = 1. / 9. * rho_i;
	mixed_term_shifted = prefactor * (u_eq_shifted0);
	mixed_term_eq = prefactor * (u_0_i);
	s_i[1] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[3] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (u_eq_shifted1);
	mixed_term_eq = prefactor * (u_1_i);
	s_i[2] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[4] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	// directions 5-8
	weighting = 1. / 36. * rho_i;
	mixed_term_shifted = prefactor * (u_eq_shifted0 + u_eq_shifted1);
	mixed_term_eq = prefactor * (u_0_i + u_1_i);
	s_i[5] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[7] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (-u_eq_shifted0 + u_eq_shifted1);
	mixed_term_eq = prefactor * (-u_0_i + u_1_i);
	s_i[6] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[8] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	// Add Source term
	f_i[0] += s_i[0];
	f_i[1] += s_i[1];
	f_i[2] += s_i[2];
	f_i[3] += s_i[3];
	f_i[4] += s_i[4];
	f_i[5] += s_i[5];
	f_i[6] += s_i[6];
	f_i[7] += s_i[7];
	f_i[8] += s_i[8];
}

} /* namespace ExternalForceFunctions */
} /* namespace natrium */

#endif /* EXTERNALFORCEFUNCTIONS_H_ */
