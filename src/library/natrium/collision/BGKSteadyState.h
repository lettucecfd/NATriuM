/**
 * @file BGKSteadyState.h
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BGKSTEADYSTATE_H_
#define BGKSTEADYSTATE_H_

#include "BGK.h"

#include <cassert>

#include "deal.II/base/index_set.h"

#include "../utilities/BasicNames.h"

#include "../utilities/Math.h"

namespace natrium {

/**
 * @short BGK model for preconditioned Navier-Stokes equation
 *        Paper by Guo et al. (2004): Preconditioned lattice-Boltzmann
 *        method for steady flows, Physical Review E 70.
 *        The equilibrium distribution functions agrees with the one
 *        of the standard BGK model, except it has the factor
 *        \f[\frac{1}{\gamma}\f] with each \f[ O(u^2) \f] term.
 *		  The hydrodynamic equations of the Steady State BGK can be
 *		  derived by the Chapman-Enskog analysis. Except for the
 *		  temporal derivative, the compressible Navier-Stokes equation
 *		  is found, but with a much better eigenvalue system. This makes
 *		  the preconditioned BGK perfect for steady flows, although not physical
 *		  for unsteady flows.
 */
class BGKSteadyState: public BGK {
private:
	double m_gamma;
public:


	/** @short Constructor
	 *  @param preconditioning_parameter (\f[ \gamma \in (0,1] \f] is the preconditioning parameter
	 *         \f[ \gamma = 1 \f] gives the standard BGK model, for  \f[ \gamma \leftarrow 0 \f],
	 *        the effective sound speed is decreased
	 *
	 */
	BGKSteadyState(double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil, double preconditioning_parameter);


	/// destructor
	virtual ~BGKSteadyState();

	/** @short function for the calculation of the equilibrium distribution in the incompressible D2Q9 model
	 *  @param i index of the direction
	 *  @param u macroscopic velocity
	 *  @param rho macroscopic density
	 *  @return value of the equilibrium distribution
	 *  @note The calculation can surely be done more efficiently by passing different arguments,
	 *        e.g. u*u or u/(c^2)
	 */
	virtual double getEquilibriumDistribution(size_t i,
			const numeric_vector& u, const double rho = 1) const;


	/**
	 * @short function for collision
	 * @short f the global vectors of discrete particle distribution functions
	 * @short densities the global vector of densities
	 * @short velocities the global vectors of velocity components [ [u_1x, u_2x, ...], [u_1y, u_2y, ...] ]
	 * @short inInitializationProcedure indicates if the collision is performed in the context of an iterative initilizatation procedure. In this case, only the macroscopic densities are recalculated, while the velocities remain unchanged. default: false
	 */
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	/**
	 * @short calculate relaxation parameter
	 * @note the preconditioning parameter gamma is only used for steady state
	 */
	static double calculateRelaxationParameter(double viscosity,
			double timeStepSize, const Stencil& stencil, double preconditioning_parameter) {
		assert(viscosity > 0.0);
		assert(timeStepSize > 0.0);
		return (viscosity) / (timeStepSize * stencil.getSpeedOfSoundSquare() * preconditioning_parameter);
	}

};


} /* namespace natrium */
#endif /* BGKSTEADYSTATE_H_ */
