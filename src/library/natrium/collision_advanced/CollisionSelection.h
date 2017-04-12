/*
 * CollisionSelection.h
 *
 *  Created on: 30.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_ADVANCED_COLLISIONSELECTION_H_
#define LIBRARY_NATRIUM_COLLISION_ADVANCED_COLLISIONSELECTION_H_
#include <vector>
#include "../stencils/Stencil.h"
#include "AuxiliaryCollisionFunctions.h"
#include "Equilibria.h"
#include "CollisionSchemes.h"
#include "CollisionOperator.h"

namespace natrium{
inline void selectCollision(
		boost::shared_ptr<SolverConfiguration>& configuration,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity, const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure) {

	switch (configuration->getStencil()) {
	case Stencil_D2Q9:
		switch (configuration->getCollisionScheme()) {
		case BGK_STANDARD:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<2,9> prams(stencil.getScaling(),viscosity,stencil.getSpeedOfSoundSquare(),delta_t);

				CollisionOperator<2,9,BGKEquilibrium,BGKCollision> BGK_BGKEQ_D2Q9;

				BGK_BGKEQ_D2Q9.collideAll(f, densities, velocities, locally_owned_dofs,
						inInitializationProcedure, prams);
			}
		}

	}
}
}



#endif /* LIBRARY_NATRIUM_COLLISION_ADVANCED_COLLISIONSELECTION_H_ */
