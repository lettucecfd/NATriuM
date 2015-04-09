/**
 * @file CollisionModel.h
 * @short Abstract class for the description of collision schemes.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COLLISIONMODEL_H_
#define COLLISIONMODEL_H_

#include <exception>

#include "boost/shared_ptr.hpp"

#include "../stencils/Stencil.h"

#include "../solver/DistributionFunctions.h"

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"

namespace natrium {

/**
 * @short Exception class for Collision
 */
class CollisionException: public NATriuMException {
private:
	std::string message;
public:
	CollisionException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	CollisionException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~CollisionException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short Abstract collision model. Required to have a common parent of all template specializations of Collision
 */
class CollisionModel {
private:
	shared_ptr<Stencil> m_stencil;
public:
	CollisionModel(const shared_ptr<Stencil> stencil) :
			m_stencil(stencil) {
	}
	;
	virtual ~CollisionModel() {
	}
	;
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const = 0;
	virtual void setTimeStep(double dt) = 0;

	const shared_ptr<Stencil>& getStencil() const {
		return m_stencil;

	}


	//////////////////////////////
	// EQUILIBRIUM DISTRIBUTION //
	//////////////////////////////

	/** @short virtual function for the calculation of the equilibrium distribution
	 *  @param i index of the direction
	 *  @param u macroscopic velocity
	 *  @param rho macroscopic density
	 *  @return value of the equilibrium distribution
	 *  @note The calculation can surely be done more efficiently by passing different arguments,
	 *        e.g. u*u or u/(c^2)
	 */
	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const = 0;


	/** @short function for the calculation of all equilibrium distributions
	 *  @param[out] feq vector of all equality distributions, must have size Q
	 *  @param[in] u macroscopic velocity
	 *  @param[in] rho macroscopic density
	 *  @note The calculation can surely be done more efficiently by passing different arguments,
	 *        e.g. u*u or u/(c^2)
	 */
	virtual void getEquilibriumDistributions(vector<double>& feq,
			const numeric_vector& u, const double rho = 1) const;
};

} /* namespace natrium */
#endif /* COLLISIONMODEL_H_ */
