/*
 * MRTEntropic.h
 *
 *  Created on: 29.03.2016
 *      Author: Dominik Wilde
 */

#ifndef MRTENTROPIC_H_
#define MRTENTROPIC_H_

#include "MRTStandard.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

class MRTEntropic: public MRT {
public:
	MRTEntropic(double relaxationParameter, double dt,
			const boost::shared_ptr<Stencil> stencil);
	virtual ~MRTEntropic();
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	/**
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD2Q9(DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;

	void collideAllD3Q19(DistributionFunctions& f, distributed_vector& densities,
				vector<distributed_vector>& velocities,
				const dealii::IndexSet& locally_owned_dofs,
				bool inInitializationProcedure) const;

	};





}


#endif /* MRTEntropic_H_ */
