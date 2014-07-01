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
class CollisionModelException: public std::exception {
private:
	std::string message;
public:
	CollisionModelException(const char *msg) :
			message(msg) {
	}
	~CollisionModelException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};


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
			DistributionFunctions& f) const = 0;


	/**
	 * @short function for collision
	 * @short f the global vectors of discrete particle distribution functions
	 * @short densities the global vector of densities
	 * @short velocities the global vectors of velocity components [ [u_1x, u_2x, ...], [u_1y, u_2y, ...] ]
	 * @short inInitializationProcedure indicates if the collision is performed in the context of an iterative initilizatation procedure. In this case, only the macroscopic densities are recalculated, while the velocities remain unchanged. default: false
	 */
	void collideAll(DistributionFunctions& f,
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
	void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const;

	/**
	 * @short calculate relaxation parameter
	 */
	static double calculateRelaxationParameter(double viscosity,
			double timeStepSize,
			const boost::shared_ptr<BoltzmannModel> boltzmannModel) {
		assert(viscosity > 0.0);
		assert(timeStepSize > 0.0);
		return (viscosity)
				/ (timeStepSize * boltzmannModel->getSpeedOfSoundSquare());
	}

	/**
	 * @short calculate the time step, so that the Courant number is 1 for the diagonal directions
	 */
	static double calculateOptimalTimeStep(double dx,
			const boost::shared_ptr<BoltzmannModel> boltzmannModel) {
		assert(dx > 0);
		return dx / boltzmannModel->getMaxParticleVelocityMagnitude();
	}

};

} /* namespace natrium */
#endif /* COLLISIONMODEL_H_ */
