/*
 * ConfigNames.h
 *
 *  Created on: 26.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_UTILITIES_CONFIGNAMES_H_
#define LIBRARY_NATRIUM_UTILITIES_CONFIGNAMES_H_

namespace natrium {

//////////////////////////////
// DECLARE ALL SELECTIONS ////
//////////////////////////////

/**
 * @short relaxation matrix for MRT
 */
enum RelaxMode {
	RELAX_FULL, // relax all higher-order moments to equilibrium
	DELLAR_RELAX_ONLY_N // relax only J to equilibrium
};

/**
 * @short moment basis for MRT
 */
enum MomentBasis {
	DELLAR_D2Q9,
	LALLEMAND_D2Q9,
	DHUMIERES_D3Q19
};


/**
 * @short Implemented streaming data types
 */
enum AdvectionSchemeName {
	SEDG, SEMI_LAGRANGIAN
};

/**
 * @short Implemented collision models
 */
enum CollisionSchemeName {
	BGK_STANDARD, // Standard BGK collision Collision for the distribution function as defined in MinLee2011
	BGK_STANDARD_TRANSFORMED, // BGK collisions with transformed distributions, as used in Palabos
	BGK_STEADY_STATE, // Steady state preconditioning by Guo et al. (2004)
	BGK_MULTIPHASE,
	BGK_INCOMPRESSIBLE, // BGK collision for incompressible Navier Stokes equations by He & Luo (1997)
	MRT_STANDARD, // Multiple Relaxation Time scheme by d'Humières (1992)
	MRT_ENTROPIC, // Optimized Multiple Relaxation Time scheme by A.Kraemer (2016)
	ENTROPIC_STABILIZED, // Optimized Multiple Relaxation Time scheme by A.Kraemer (2016)
	KBC_STANDARD, // Multiple Relaxation Time scheme with autonomously adaptive parameters by Karlin et al. (2014)
	KBC_CENTRAL, // // Multiple Relaxation Time scheme with autonomously adaptive parameters by Karlin et al. (2014), central moments are used
	BGK_MULTI_AM4, // Multistep BGK model according to Krämer (2017, Uni Siegen, Anhang A)
	BGK_MULTI_BDF2,
	BGK_REGULARIZED // Regularized model by Latt and Chopard (2005)
};

enum EquilibriumSchemeName {
	BGK_EQUILIBRIUM
};


// StencilType defined in Stencil.h

/**
 * @short Implemented time integrators
 * @note other refers to dealii integrators which are accessed through a wrapper class
 */
enum TimeIntegratorName {
	RUNGE_KUTTA_5STAGE, THETA_METHOD, EXPONENTIAL, OTHER
};

enum DealIntegratorName {
	FORWARD_EULER,
	RK_THIRD_ORDER,
	RK_CLASSIC_FOURTH_ORDER,
	BACKWARD_EULER,
	IMPLICIT_MIDPOINT,
	CRANK_NICOLSON,
	SDIRK_TWO_STAGES,
	HEUN_EULER,
	BOGACKI_SHAMPINE,
	DOPRI,
	FEHLBERG,
	CASH_KARP,
	NONE
};

enum DealSolverName {
	BICGSTAB, CG, FGMRES, GMRES, MINRES, QMRS, RELAXATION, RICHARDSON
};

/**
 * @short the numerical flux used to calculate the advection operator
 */
enum FluxTypeName {
	LAX_FRIEDRICHS, CENTRAL
};

/**
 * @short the initialization procedure for the distribution functions
 */
enum InitializationSchemeName {
	EQUILIBRIUM, // Distribute with equilibrium functions
	ITERATIVE // Distribute with iterative procedure; enforces consistent initial conditions
};

enum PseudopotentialType {
	SHAN_CHEN, SUKOP, CARNAHAN_STARLING
};

enum ForceType {
	NO_FORCING, 	// No external force
	SHIFTING_VELOCITY,  // Shifting velocity method
	EXACT_DIFFERENCE,   // Exact difference method by Kuppershtokh
	GUO
};

enum FilteringSchemeName {
	EXPONENTIAL_FILTER, NEW_FILTER
};

enum SupportPointsName {
	GAUSS_LOBATTO_POINTS, GAUSS_CHEBYSHEV_POINTS, GAUSS_LOBATTO_CHEBYSHEV_POINTS, EQUIDISTANT_POINTS
};

enum QuadratureName {
	QGAUSS_LOBATTO, QGAUSS
};

enum RegularizationSchemeName {
	NO_REGULARIZATION, PSEUDO_ENTROPY_MAXIMIZATION, ZERO_HIGH_ORDER_MOMENTS, ENTROPY_MAXIMIZATION, PSEUDO_ENTROPY_MAXIMIZATION_WITH_E
};

}


#endif /* LIBRARY_NATRIUM_UTILITIES_CONFIGNAMES_H_ */
