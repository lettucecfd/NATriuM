/**
 * @file CollisionModel.h
 * @short Abstract class for the description of collision schemes.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COLLISIONMODEL_H_
#define COLLISIONMODEL_H_

#include "boost/shared_ptr.hpp"

#include "../boltzmannmodels/BoltzmannModel.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Abstract class for the description of collision schemes.
 */
class CollisionModel {

protected:
	/// Boltzmann model (e.g. D2Q9Incompressible)
	const boost::shared_ptr<BoltzmannModel> m_boltzmannModel;

	/// D (dimension)
	size_t m_d;

	/// Q (number of directions)
	double m_q;

public:

	/**
	 * @short constructor
	 * @param[in] viscosity the kinematic viscosity
	 * @param[in] timeStepSize the size of the time step
	 * @param[in] boltzmannModel the discrete velocity stencil
	 */
	CollisionModel(boost::shared_ptr<BoltzmannModel> boltzmannModel);

	/// destructor
	virtual ~CollisionModel();

	/**
	 * @short virtual function for collision
	 * @param[in/out] distributions the particle distribution functions
	 */
	virtual void collideSinglePoint(vector<double>& distributions) const = 0;

	/**
	 * @short virtual function for collision
	 * @param[in] doF the doF index for which collision is done
	 * @param[in] feq the vector of local equilibrium distributions
	 * @param[in] f the vector of global distribution functions
	 */
	virtual void collideSingleDoF(size_t doF, const vector<double>& feq,
			vector<distributed_vector>& f) const = 0;

	/**
	 * @short function for collision
	 * @short f the global vectors of discrete particle distribution functions
	 * @short densities the global vector of densities
	 * @short velocities the global vectors of velocity components [ [u_1x, u_2x, ...], [u_1y, u_2y, ...] ]
	 */
	void collideAll(vector<distributed_vector>& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities) const;

};

} /* namespace natrium */
#endif /* COLLISIONMODEL_H_ */
