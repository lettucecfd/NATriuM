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
			const boost::shared_ptr<Stencil> stencil);
	virtual ~MRT();

	// virtual void collideSinglePoint(vector<double>& distributions) const;

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho) const;


	void setTimeStep(double dt) {
		assert(dt > 0);
		double tau_times_dt = m_dt * m_relaxationParameter;
		m_dt = dt;
		CollisionModel::setTimeStep(dt);
		m_relaxationParameter = tau_times_dt / dt;
		m_prefactor = -1. / (m_relaxationParameter + 0.5);
	}

//	double getTimeStep() const
//	{
//		return m_dt;
//	}

	void setRelaxationParameter(double tau, double dt) {
		assert(tau > 0);
		m_relaxationParameter = tau;
		//cout << m_relaxationParameter << endl;
		m_prefactor = -1. / (tau + 0.5);
		m_dt = dt;
	}

	double getTime() const
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

