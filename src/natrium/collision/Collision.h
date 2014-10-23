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

#include "../boltzmannmodels/BoltzmannModel.h"

#include "../solver/DistributionFunctions.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short Exception class for Collision model
 */
class CollisionException: public std::exception {
private:
	std::string message;
public:
	CollisionException(const char *msg) :
			message(msg) {
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
class CollisionModel{
public:
	CollisionModel(){};
	virtual ~CollisionModel(){};
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities, const vector<bool>& isBoundary,
			bool inInitializationProcedure = false) const = 0;
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const = 0;
};

/** @short Abstract class for the description of collision schemes.
 */
template <class BoltzmannType, class CollisionType>
class Collision: public CollisionModel {

protected:
	/// Boltzmann model (e.g. D2Q9Incompressible)
	const BoltzmannType m_boltzmannType;

	const CollisionType m_collisionType;

	/// D (dimension)
	size_t m_d;

	/// Q (number of directions)
	double m_q;

public:

	/**
	 * @short constructor, important: Copy arguments instead of reference
	 * @param[in] viscosity the kinematic viscosity
	 * @param[in] timeStepSize the size of the time step
	 * @param[in] boltzmannModel the discrete velocity stencil
	 */
	Collision(BoltzmannType boltzmannType, CollisionType collisionType) ;

	/// destructor
	virtual ~Collision();

	/**
	 * @short only for testing purposes; do not use in code, because this function is very inefficient
	 */
	virtual void collideSinglePoint(vector<double>& distributions) const {

		// assert
		assert(distributions.size() == m_q);

		// create DistributionFunctions
		DistributionFunctions f;
		f.reinit(m_q, 1);
		for (size_t i = 0; i < m_q; i++){
			f.at(i)(0) = distributions.at(i);
		}


		// calculate macroscopic entities
		double rho = m_boltzmannType.calculateDensity(distributions);
		numeric_vector u(m_d);
		m_boltzmannType.calculateVelocity(distributions, rho, u);

		// calculate equilibrium distribution (feq)
		vector<double> feq(m_q);
		m_boltzmannType.getEquilibriumDistributions(feq, u, rho);

		// update distribution
		m_collisionType.collideSingleDoF(0, feq, f);
		for (size_t i = 0; i < m_q; i++) {
			distributions.at(i) += m_collisionType.getPrefactor()
					* (distributions.at(i) - feq.at(i));
		}

	} //collide

	/**
	 * @short function for collision
	 * @short f the global vectors of discrete particle distribution functions
	 * @short densities the global vector of densities
	 * @short velocities the global vectors of velocity components [ [u_1x, u_2x, ...], [u_1y, u_2y, ...] ]
	 * @short inInitializationProcedure indicates if the collision is performed in the context of an iterative initilizatation procedure. In this case, only the macroscopic densities are recalculated, while the velocities remain unchanged. default: false
	 */
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities, const vector<bool>& isBoundary,
			bool inInitializationProcedure = false) const;

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


	/**
	 * @short calculate the time step, so that the Courant number is 1 for the diagonal directions
	 */
	static double calculateOptimalTimeStep(double dx,
			const BoltzmannType& boltzmannType) {
		assert(dx > 0);
		return dx / boltzmannType.getMaxParticleVelocityMagnitude();
	}

};

} /* namespace natrium */
#endif /* COLLISIONMODEL_H_ */
