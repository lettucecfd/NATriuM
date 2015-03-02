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

#include <cassert>

#include "../utilities/Math.h"

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
			const double rho = 1) const{
		// TODO efficient implementation of equilibrium distribution (for whole vector of equilibrium distributions)
		// TODO efficient implementation of scalar product between vectors and 'unit vectors' like (-1,-1); case i=...
		// TODO save 1./(speedOfSoundSquare) in private variable
		// TODO /= 2 by bitshift
		// TODO (..)-Term in (9.19) -> calculate only once and save in private variable
		//vector_expression<double>::inner_prod_prec(u,u);
		assert(i < Q);
		assert(rho > 0);
		assert(i >= 0);
		assert(u.size() == D);
		assert(u(0) < 1000000000000000.);
		assert(u(1) < 1000000000000000.);

		double prefactor = getWeight(i) * rho;
		double uSquareTerm  = - Math::scalar_product(u, u)/(2*m_speedOfSoundSquare);
		if (0 == i){
			return prefactor * (1 + uSquareTerm);
		}
		double mixedTerm = Math::scalar_product(u,getDirection(i))/m_speedOfSoundSquare;
		return prefactor * (1 + mixedTerm * (1 + 0.5*(mixedTerm)) + uSquareTerm);

	}/// getEquilibriumDistribution



	// TODO getEquilibriumDistributions implementation -> more efficient

};

} /* namespace natrium */
#endif /* D2Q9INCOMPRESSIBLE_H_ */
