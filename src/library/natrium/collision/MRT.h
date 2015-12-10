/*
 * MRT.h
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#ifndef MRT_H_
#define MRT_H_

#include "CollisionModel.h"
#include "../utilities/BasicNames.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

class MRT: public CollisionModel {
private:

	/// relaxation time Tau
	double m_relaxationParameter;

	/// prefactor of the collision
	double m_prefactor;

	// time step size
	double m_dt;

public:
	MRT(double relaxationParameter, double dt,
			const shared_ptr<Stencil> stencil);
	virtual ~MRT();

	// virtual void collideSinglePoint(vector<double>& distributions) const;

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho) const;

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const = 0;

	static double calculateRelaxationParameter(double viscosity,
			double timeStepSize, const Stencil& stencil,
			double preconditioning_parameter = 1.0) {
		assert(viscosity > 0.0);
		assert(timeStepSize > 0.0);
		return (viscosity) / (timeStepSize * stencil.getSpeedOfSoundSquare());
	}

	void setTimeStep(double dt) {
		assert(dt > 0);
		double tau_times_dt = m_dt * m_relaxationParameter;
		m_dt = dt;
		m_relaxationParameter = tau_times_dt / dt;
		m_prefactor = -1. / (m_relaxationParameter + 0.5);
	}

	double getTimeStep() const
	{
		return m_dt;
	}

	size_t getQ() const {
		return getStencil()->getQ();
	}

	double getPrefactor() const {
		return m_prefactor;
	}

	double getRelaxationParameter() const {
		return m_relaxationParameter;
	}


};
}
// namespace natrium

#endif /* MRT_H_ */

