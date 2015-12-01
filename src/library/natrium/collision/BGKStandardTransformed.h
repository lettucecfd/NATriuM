/**
 * @file BGKStandardTransformed.h
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BGKSTANDARDTRANSFORMED_H_
#define BGKSTANDARDTRANSFORMED_H_

#include "BGK.h"

#include <cassert>

#include "deal.II/base/index_set.h"

#include "../utilities/BasicNames.h"

#include "../utilities/Math.h"

namespace natrium {

/**
 * @short Simple BGK model with transformed distributions in order to lower round-off errors
 *        The scheme is described in "How to improve the accuracy of Lattice Boltzmann calculations"
 *        by Bastien Chopard (May 2008), available on LBMethod.org
 */
class BGKStandardTransformed: public BGK {
private:
	const double m_rho0;
public:


	/// constructor
	BGKStandardTransformed(double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil);


	/// destructor
	virtual ~BGKStandardTransformed();

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

	double getRho0() const {
		return m_rho0;
	}

	/**
	 * @short calculate macroscopic density
	 * @param[in] distributions particle distribution functions at a given point
	 * @return macroscopic density (sum of all distributions)
	 */
	virtual double calculateDensity(const vector<double>& distributions) const {

		// calculate macroscopic density (rho)
		double rho = 0.0;
		for (size_t i = 0; i < getStencil()->getQ(); i++) {
			rho += distributions.at(i);
		}
		return rho + m_rho0;

	}


};


} /* namespace natrium */
#endif /* BGKSTANDARDTRANSFORMED_H_ */
