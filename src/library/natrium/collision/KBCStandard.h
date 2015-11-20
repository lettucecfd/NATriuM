/*
 * KBCStandard.h
 *
 *  Created on: 17.11.2015
 *      Author: dominik
 */

#ifndef KBCSTANDARD_H_
#define KBCSTANDARD_H_

#include "MRT.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

class KBCStandard: public MRT {
public:
	KBCStandard(double relaxationParameter, double dt,
			const shared_ptr<Stencil> stencil);
	virtual ~KBCStandard();

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure) const;

};
} /* namespace natrium */

#endif /* KBCSTANDARD_H_ */
