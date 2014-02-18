/**
 * @file D2Q9Incompressible.h
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef D2Q9INCOMPRESSIBLEMODEL_H_
#define D2Q9INCOMPRESSIBLEMODEL_H_

#include "D2Q9Model.h"

#include "../utilities/BasicNames.h"


namespace natrium {

/**
 * @short D2Q9 model description for incompressible flow.
 */
class D2Q9IncompressibleModel: public D2Q9Model {
public:


	/// constructor
	D2Q9IncompressibleModel(double scaling = 1);


	/// destructor
	virtual ~D2Q9IncompressibleModel();


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


	// TODO getEquilibriumDistributions implementation -> more efficient

};

} /* namespace natrium */
#endif /* D2Q9INCOMPRESSIBLE_H_ */
