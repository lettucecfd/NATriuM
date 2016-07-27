/*
 * BGKMultistep.h
 *
 *  Created on: 20.06.2016
 *      Author: dominik
 */

#ifndef BGKMULTISTEP_H_
#define BGKMULTISTEP_H_

#include "MultistepCollisionData.h"
#include "BGK.h"
#include "ExternalForceFunctions.h"
#include <cassert>

#include "deal.II/base/index_set.h"

#include "../utilities/BasicNames.h"

#include "../utilities/Math.h"

namespace natrium {

class BGKMultistep: public BGK, public MultistepCollisionData  {
public:
	enum MultistepModelName
	{
		ADAMSMOULTON4,
		BDF2
	};



	BGKMultistep(double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil, int model);

	MultistepModelName m_model;

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const;

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	/**
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD2Q9(DistributionFunctions& f,
			distributed_vector& densities, vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;

	void collideAllD3Q19(DistributionFunctions& f,
			distributed_vector& densities, vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;



	virtual ~BGKMultistep();
};

} /* namespace natrium */

#endif /* BGKMULTISTEP_H_ */
