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

	/// relaxation parameter
	double m_relaxationParameter;

	/// prefactor of the collision (- 1/(tau + 0.5))
	double m_prefactor;

	// time step size
	double m_dt;

public:
	MRT(double relaxationParameter,double dt, const shared_ptr<Stencil> stencil);
	virtual ~MRT();

	virtual void collideSinglePoint(vector<double>& distributions) const;

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const;

	void collideSingleDoF(size_t doF, const vector<double>& feq,
			DistributionFunctions& f) const {
		for (size_t j = 0; j < getStencil()->getQ(); j++) {
			f.at(j)(doF) += m_prefactor * (f.at(j)(doF) - feq.at(j));
		}
	}

	void setTimeStep(double dt) {
		assert(dt > 0);
		m_dt = dt;
	}

	size_t getQ() const {
		return getStencil()->getQ();
	}

	double getPrefactor() const {
		return m_prefactor;
	}


};
}
// namespace natrium

#endif /* MRT_H_ */

