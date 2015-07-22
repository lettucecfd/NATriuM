/*
 * BGKIncompressible.h
 *
 *  Created on: 19.07.2015
 *      Author: dominik
 */

#ifndef BGKINCOMPRESSIBLE_H_
#define BGKINCOMPRESSIBLE_H_

#include "../collision/BGK.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"
namespace natrium {

class BGKIncompressible: public BGK {
public:
	BGKIncompressible(double relaxationParameter, double dt,
			const shared_ptr<Stencil> stencil);
	virtual ~BGKIncompressible();

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const;

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false) const;

};

} /* namespace natrium */

#endif /* BGKINCOMPRESSIBLE_H_ */
