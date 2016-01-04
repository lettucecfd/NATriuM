/**
 * @file BGKIncompressible.h
 * @short D2Q9 model description for improved incompressible flow by Guo et al.(2000)s.
 * @date 19.07.2015
 * @author Dominik Wilde, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BGKINCOMPRESSIBLE_H_
#define BGKINCOMPRESSIBLE_H_

#include "BGK.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"
namespace natrium {

class BGKIncompressible: public BGK {
private:
	double m_lambda = 1.0 / 3.0;
	double m_gamma = 1.0 / 12.0;
public:


	BGKIncompressible(double relaxationParameter, double dt,
			const boost::shared_ptr<Stencil> stencil);
	virtual ~BGKIncompressible();

	/** @short function for the calculation of the equilibrium distribution in the incompressible D2Q9 model
	 *  @param i index of the direction
	 *  @param u macroscopic velocity
	 *  @param rho macroscopic density
	 *  @return value of the equilibrium distribution
	 *  @note The calculation can surely be done more efficiently by passing different arguments,
	 *        e.g. u*u or u/(c^2)
	 */
	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const;

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
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD2Q9(DistributionFunctions& f,
			distributed_vector& densities, vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;


};

} /* namespace natrium */

#endif /* BGKINCOMPRESSIBLE_H_ */
