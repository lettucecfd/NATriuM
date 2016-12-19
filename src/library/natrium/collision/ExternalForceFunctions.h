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

/**
 * @short Calculate source term for exact difference forcing (Kupershtokh et. al.)
 * @note See Phys.Rev. E 84, 046710 (2011), Scheme V
 */
inline void applyExactDifferenceForcingD2Q9(double f_i[], double force_x,
		double force_y, double u_0_i, double u_1_i, double rho_i, double dt,
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
	double uSquareTerm_eq_shifted = -0.5 * prefactor
			* scalar_product_eq_shifted;

	double scalar_product_eq = u_0_i * u_0_i + u_1_i * u_1_i;
	double u_square_term_eq = -0.5 * prefactor * scalar_product_eq;

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

inline void applyGuoForcingD2Q9(double f_i[], double force_x, double force_y,
		double u_0_i, double u_1_i, double omega, double prefactor, double dt) {

	double s_i[9];

	double eia_ua;
	double eig_fg;
	double ug_fg = u_0_i * force_x + u_1_i * force_y;

	double wi_one_min_omega_half;
	double one_min_omega_half = 1 - 0.5 * omega;

	//f0
	eia_ua = 0;
	eig_fg = 0;
	wi_one_min_omega_half = 4. / 9. * one_min_omega_half;
	s_i[0] = wi_one_min_omega_half * (prefactor * (eig_fg - ug_fg));
	//f1
	wi_one_min_omega_half = 1. / 9. * one_min_omega_half;
	eia_ua = u_0_i;
	eig_fg = force_x;
	s_i[1] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f2
	eia_ua = u_1_i;
	eig_fg = force_y;
	s_i[2] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f3
	eia_ua = -u_0_i;
	eig_fg = -force_x;
	s_i[3] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f4
	eia_ua = -u_1_i;
	eig_fg = -force_y;
	s_i[4] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f5
	wi_one_min_omega_half = 1. / 36. * one_min_omega_half;
	eia_ua = u_0_i + u_1_i;
	eig_fg = force_x + force_y;
	s_i[5] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f6
	eia_ua = -u_0_i + u_1_i;
	eig_fg = -force_x + force_y;
	s_i[6] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f7
	eia_ua = -u_0_i - u_1_i;
	eig_fg = -force_x - force_y;
	s_i[7] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f8
	eia_ua = u_0_i - u_1_i;
	eig_fg = force_x - force_y;
	s_i[8] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));

	// Add Source term
	f_i[0] += dt * s_i[0];
	f_i[1] += dt * s_i[1];
	f_i[2] += dt * s_i[2];
	f_i[3] += dt * s_i[3];
	f_i[4] += dt * s_i[4];
	f_i[5] += dt * s_i[5];
	f_i[6] += dt * s_i[6];
	f_i[7] += dt * s_i[7];
	f_i[8] += dt * s_i[8];
}

/**
 * @short Calculate source term for exact difference forcing (Kupershtokh et. al.)
 * @note See Phys.Rev. E 84, 046710 (2011), Scheme V
 */
inline void applyExactDifferenceForcingD3Q19(double f_i[], double force_x,
		double force_y, double force_z, double u_0_i, double u_1_i,
		double u_2_i, double rho_i, double dt, double prefactor) {
	double mixed_term_shifted;
	double mixed_term_eq;
	double s_i[19];

	// Calculate velocity shift
	double deltaU_0 = force_x * dt / rho_i;
	double deltaU_1 = force_y * dt / rho_i;
	double deltaU_2 = force_z * dt / rho_i;

	double u_eq_shifted0 = u_0_i + deltaU_0;
	double u_eq_shifted1 = u_1_i + deltaU_1;
	double u_eq_shifted2 = u_2_i + deltaU_2;

	double scalar_product_eq_shifted = u_eq_shifted0 * u_eq_shifted0
			+ u_eq_shifted1 * u_eq_shifted1 + u_eq_shifted2 * u_eq_shifted2;
	double uSquareTerm_eq_shifted = -0.5 * prefactor
			* scalar_product_eq_shifted;

	double scalar_product_eq = u_0_i * u_0_i + u_1_i * u_1_i + u_2_i * u_2_i;
	double u_square_term_eq = -0.5 * prefactor * scalar_product_eq;

	// Calculate source term
	// direction 0
	double weighting = 1. / 3. * rho_i;
	s_i[0] = weighting * (1 + uSquareTerm_eq_shifted)
			- weighting * (1 + u_square_term_eq);

	// directions 1-6
	weighting = 1. / 18. * rho_i;
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
	s_i[6] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[5] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (u_eq_shifted2);
	mixed_term_eq = prefactor * (u_2_i);
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
	// directions 7-18
	weighting = 1. / 36. * rho_i;
	mixed_term_shifted = prefactor * (u_eq_shifted0 + u_eq_shifted1);
	mixed_term_eq = prefactor * (u_0_i + u_1_i);
	s_i[12] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[14] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (u_eq_shifted0 - u_eq_shifted1);
	mixed_term_eq = prefactor * (u_0_i - u_1_i);
	s_i[11] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[13] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (u_eq_shifted0 + u_eq_shifted2);
	mixed_term_eq = prefactor * (u_0_i + u_2_i);
	s_i[7] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[9] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (u_eq_shifted0 - u_eq_shifted2);
	mixed_term_eq = prefactor * (u_0_i - u_2_i);
	s_i[10] = weighting
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
	mixed_term_shifted = prefactor * (u_eq_shifted1 + u_eq_shifted2);
	mixed_term_eq = prefactor * (u_1_i + u_2_i);
	s_i[16] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[18] = weighting
			* (1 - mixed_term_shifted * (1 - 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 - mixed_term_eq * (1 - 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	mixed_term_shifted = prefactor * (u_eq_shifted1 - u_eq_shifted2);
	mixed_term_eq = prefactor * (u_1_i - u_2_i);
	s_i[17] = weighting
			* (1 + mixed_term_shifted * (1 + 0.5 * mixed_term_shifted)
					+ uSquareTerm_eq_shifted)
			- weighting
					* (1 + mixed_term_eq * (1 + 0.5 * mixed_term_eq)
							+ u_square_term_eq);
	s_i[15] = weighting
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
	f_i[9] += s_i[9];
	f_i[10] += s_i[10];
	f_i[11] += s_i[11];
	f_i[12] += s_i[12];
	f_i[13] += s_i[13];
	f_i[14] += s_i[14];
	f_i[15] += s_i[15];
	f_i[16] += s_i[16];
	f_i[17] += s_i[17];
	f_i[18] += s_i[18];
}

inline void applyGuoForcingD3Q19(double f_i[], double force_x, double force_y,
		double force_z, double u_0_i, double u_1_i, double u_2_i, double omega,
		double prefactor, double dt) {

	double s_i[19];

	double eia_ua;
	double eig_fg;
	double ug_fg = u_0_i * force_x + u_1_i * force_y + u_2_i * force_z;

	double wi_one_min_omega_half;
	double one_min_omega_half = 1 - 0.5 * omega;

	//f0
	eia_ua = 0;
	eig_fg = 0;
	wi_one_min_omega_half = 1. / 3. * one_min_omega_half;
	s_i[0] = wi_one_min_omega_half * (prefactor * (eig_fg - ug_fg));
	//f1
	wi_one_min_omega_half = 1. / 18. * one_min_omega_half;
	eia_ua = u_0_i;
	eig_fg = force_x;
	s_i[1] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f2
	eia_ua = u_2_i;
	eig_fg = force_z;
	s_i[2] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f3
	eia_ua = -u_0_i;
	eig_fg = -force_x;
	s_i[3] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f4
	eia_ua = -u_2_i;
	eig_fg = -force_z;
	s_i[4] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f5
	eia_ua = -u_1_i;
	eig_fg = -force_y;
	s_i[5] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f6
	eia_ua = u_1_i;
	eig_fg = force_y;
	s_i[6] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f7
	wi_one_min_omega_half = 1. / 36. * one_min_omega_half;
	eia_ua = u_0_i + u_2_i;
	eig_fg = force_x + force_z;
	s_i[7] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f8
	eia_ua = -u_0_i + u_2_i;
	eig_fg = -force_x + force_z;
	s_i[8] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f9
	eia_ua = -u_0_i - u_2_i;
	eig_fg = -force_x - force_z;
	s_i[9] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f10
	eia_ua = u_0_i - u_2_i;
	eig_fg = force_x - force_z;
	s_i[10] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f11
	eia_ua = u_0_i - u_1_i;
	eig_fg = force_x - force_y;
	s_i[11] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f12
	eia_ua = u_0_i + u_1_i;
	eig_fg = force_x + force_y;
	s_i[12] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f13
	eia_ua = -u_0_i + u_1_i;
	eig_fg = -force_x + force_y;
	s_i[13] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f14
	eia_ua = -u_0_i - u_1_i;
	eig_fg = -force_x - force_y;
	s_i[14] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f15
	eia_ua = -u_1_i + u_2_i;
	eig_fg = -force_y + force_z;
	s_i[15] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f16
	eia_ua = u_1_i + u_2_i;
	eig_fg = force_y + force_z;
	s_i[16] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f17
	eia_ua = u_1_i - u_2_i;
	eig_fg = force_y - force_z;
	s_i[17] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));
	//f18
	eia_ua = -u_1_i - u_2_i;
	eig_fg = -force_y - force_z;
	s_i[18] = wi_one_min_omega_half
			* (prefactor * (eig_fg - ug_fg)
					+ prefactor * prefactor * (eia_ua * eig_fg));

	// Add Source term
	f_i[0] += dt * s_i[0];
	f_i[1] += dt * s_i[1];
	f_i[2] += dt * s_i[2];
	f_i[3] += dt * s_i[3];
	f_i[4] += dt * s_i[4];
	f_i[5] += dt * s_i[5];
	f_i[6] += dt * s_i[6];
	f_i[7] += dt * s_i[7];
	f_i[8] += dt * s_i[8];
	f_i[9] += dt * s_i[9];
	f_i[10] += dt * s_i[10];
	f_i[11] += dt * s_i[11];
	f_i[12] += dt * s_i[12];
	f_i[13] += dt * s_i[13];
	f_i[14] += dt * s_i[14];
	f_i[15] += dt * s_i[15];
	f_i[16] += dt * s_i[16];
	f_i[17] += dt * s_i[17];
	f_i[18] += dt * s_i[18];
}

} /* namespace ExternalForceFunctions */
} /* namespace natrium */

#endif /* EXTERNALFORCEFUNCTIONS_H_ */
