/**
 * @file CollisionModel.h
 * @short Abstract class for the description of collision schemes.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COLLISIONMODEL_H_
#define COLLISIONMODEL_H_

#include <exception>

#include "deal.II/base/index_set.h"

#include "../stencils/Stencil.h"

#include "../solver/DistributionFunctions.h"
#include "../solver/SolverConfiguration.h"

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"

#include "../problemdescription/ConstantExternalForce.h"


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
	boost::shared_ptr<Stencil> m_stencil;
	ForceType m_forceType;
	double m_forceX;
	double m_forceY;
	double m_forceZ;
public:
	CollisionModel(const boost::shared_ptr<Stencil> stencil, ForceType force_type = NO_FORCING) :
			m_stencil(stencil),
			m_forceType(force_type),
			m_forceX(0),
			m_forceY(0),
			m_forceZ(0){
	}
	;
	virtual ~CollisionModel() {
	}
	;
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const = 0;

	virtual void setTimeStep(double dt) = 0;

	const boost::shared_ptr<Stencil>& getStencil() const {
		return m_stencil;

	}

	////////////////////////////////////
	// CALCULATE MACROSCOPIC ENTITIES //
	////////////////////////////////////

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
		return rho;

	}

	/**
	 * @short calculate macroscopic velocity
	 * @param[in] distributions particle distribution functions at a given point
	 * @return macroscopic velocity
	 */
	virtual numeric_vector calculateVelocity(
			const vector<double>& distributions) const {

		numeric_vector u(getStencil()->getD());
		for (size_t i = 0; i < getStencil()->getQ(); i++) {
			// TODO efficient calculation of scalar*directions?
			Math::add_vector(u,
					Math::scalar_vector(distributions.at(i),
							getStencil()->getDirection(i)));
		}
		Math::scale_vector(1. / this->calculateDensity(distributions), u);
		return u;
	}

	/**
	 * @short calculate macroscopic velocity; saves the double calculation of the density
	 * @note more efficient
	 * @param[in] distributions particle distribution functions at a given point
	 * @param[in] rho macroscopic density
	 * @param[out] u macroscopic velocity
	 */
	virtual void calculateVelocity(const vector<double>& distributions,
			const double rho, numeric_vector& u) const {

		// assert
		assert(u.size() == m_stencil->getD());
		assert(u(0) == 0.0);
		assert(u(m_stencil->getD() - 1) == 0.0);

		for (size_t i = 0; i < m_stencil->getQ(); i++) {
			// TODO efficient calculation of scalar*directions?
			Math::add_vector(u,
					Math::scalar_vector(distributions.at(i),
							m_stencil->getDirection(i)));
		}
		Math::scale_vector(1. / rho, u);
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

	void setExternalForce(const ConstantExternalForce<2>& external_force){
		m_forceX = external_force.getForce()[0];
		m_forceY = external_force.getForce()[1];
	}
	void setExternalForce(const ConstantExternalForce<3>& external_force){
		m_forceX = external_force.getForce()[0];
		m_forceY = external_force.getForce()[1];
		m_forceZ = external_force.getForce()[2];

	}

	ForceType getForceType() const {
		return m_forceType;
	}


	double getForceX() const {
		return m_forceX;
	}

	double getForceY() const {
		return m_forceY;
	}

	double getForceZ() const {
		return m_forceZ;
	}

	void setForceType(ForceType forceType) {
		m_forceType = forceType;
	}
};

} /* namespace natrium */
#endif /* COLLISIONMODEL_H_ */
