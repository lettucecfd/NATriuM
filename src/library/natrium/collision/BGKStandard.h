/**
 * @file BGKStandard.h
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BGKSTANDARD_H_
#define BGKSTANDARD_H_

#include "BGK.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

/**
 * @short Simple BGK model
 */
class BGKStandard: public BGK {
public:


	/// constructor
	BGKStandard(double relaxationParameter, double dt, const shared_ptr<Stencil> stencil);


	/// destructor
	virtual ~BGKStandard();

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
			bool inInitializationProcedure = false) const;


};


} /* namespace natrium */
#endif /* BGKSTANDARD_H_ */
