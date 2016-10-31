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
			const boost::shared_ptr<Stencil> stencil);
	virtual ~MRTStandard();

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

	/**
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD2Q9(DistributionFunctions& f,
			distributed_vector& densities, vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;




};

} //namespace natrium

#endif /* MRTSTANDARD_H_ */

