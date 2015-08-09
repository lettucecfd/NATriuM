/*
 * MRTStandard.h
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#ifndef MRTSTANDARD_H_
#define MRTSTANDARD_H_

#include "MRT.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"
namespace natrium {

class MRTStandard: public MRT {
public:
	MRTStandard(double relaxationParameter, double dt,
			const shared_ptr<Stencil> stencil);
	virtual ~MRTStandard();

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const;

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho) const;

	double getMomentEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho) const;

};

} //namespace natrium

#endif /* MRTSTANDARD_H_ */

