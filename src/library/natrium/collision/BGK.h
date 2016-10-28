/**
 * @file BGK.h
 * @short Description of the BGK model for the transformed particle distributions,
 *        as described in Global data which is used by Min and Lee (2011): A spectral-element discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BGK_H_
#define BGK_H_

#include "CollisionModel.h"

#include "deal.II/base/index_set.h"

#include "../utilities/BasicNames.h"

#include "../solver/DistributionFunctions.h"

namespace natrium {

/** @short Description of the BGK model for the transformed particle distributions,
 *        as described in Global data which is used by Min and Lee (2011): A spectral-element discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 */
class BGK: public CollisionModel {
private:

	/// relaxation parameter
	double m_relaxationParameter;

	/// prefactor of the collision (- 1/(tau + 0.5))
	double m_prefactor;

	// time step size
	double m_dt;
public:

	/**
	 * @short constructor
	 * @param[in] relaxationParameter relaxation parameter tau
	 */
	BGK(double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil);

	/// destructor
	virtual ~BGK();

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

	/* @short virtual function for collision
	 * @param[in] doF the doF index for which collision is done
	 * @param[in] feq the vector of local equilibrium distributions
	 * @param[in] f the vector of global distribution functions
	 */
	void collideSingleDoF(size_t doF, const vector<double>& feq,
			DistributionFunctions& f) const {
		for (size_t j = 0; j < getStencil()->getQ(); j++) {
			f.at(j)(doF) += m_prefactor * (f.at(j)(doF) - feq.at(j));
		}
	}

	/**
	 * @short only for testing purposes; do not use in code, because this function is very inefficient
	 */
	virtual void collideSinglePoint(vector<double>& distributions) const;

	/**
	 * @short calculate relaxation parameter
	 * @note preconditioning_parameter is only used for steady state
	 */
	// TODO remove



	// GETTER AND SETTER
	// TODO remove
	void setRelaxationParameter(double tau, double dt) {
		assert(tau > 0);
		m_relaxationParameter = tau;
		m_prefactor = -1. / (tau + 0.5);
		m_dt = dt;
	}

	// TODO remove
	double getPrefactor() const {
		return m_prefactor;
	}

	/**
	 * @short changes time step and relaxation time at constant velocity and speed of sound
	 */
	// TODO remove
	void setTimeStep(double dt) {
		assert(dt > 0);
		double tau_times_dt = m_dt * m_relaxationParameter;
		m_dt = dt;
		CollisionModel::setTimeStep(dt);
		m_relaxationParameter = tau_times_dt / dt;
		m_prefactor = -1. / (m_relaxationParameter + 0.5);
	}

	size_t getQ() const {
		return getStencil()->getQ();
	}

	// TODO remove
	double getRelaxationParameter() const {
		return m_relaxationParameter;
	}

	// TODO remove
	double getDt() const {
		return m_dt;
	}

	// TODO remove
	double getLambda() const {
		return m_relaxationParameter*getStencil()->getSpeedOfSoundSquare()*m_dt;
	}
};

} /* namespace natrium */
#endif /* BGK_H_ */
