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

namespace natrium {
inline void selectCollision(
		boost::shared_ptr<SolverConfiguration>& configuration,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure) {

	bool hasCollided = 0;

	switch (configuration->getStencil()) {
	case Stencil_D2Q9:
		switch (configuration->getCollisionScheme()) {
		case BGK_STANDARD:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<2, 9> prams(stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);

				BGKCollision<2, 9, BGKEquilibrium>::uniqueData data(prams);

				CollisionOperator<2, 9, BGKEquilibrium, BGKCollision> BGK_BGKEQ_D2Q9;

				BGK_BGKEQ_D2Q9.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, prams,
						data);
				hasCollided = 1;
				break;
			}
			break;
		case BGK_REGULARIZED:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<2, 9> prams(stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);

				Regularized<2, 9, BGKEquilibrium>::uniqueData data(prams);

				CollisionOperator<2, 9, BGKEquilibrium, Regularized> REG_BGKEQ_D2Q9;

				REG_BGKEQ_D2Q9.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, prams,
						data);
				hasCollided = 1;
				break;
			}
			break;
		}
		break;

	case Stencil_D3Q19:
		switch (configuration->getCollisionScheme()) {
		case BGK_STANDARD:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<3, 19> prams(stencil.getScaling(),
						viscosity, stencil, stencil.getSpeedOfSoundSquare(),
						delta_t);

				BGKCollision<3, 19, BGKEquilibrium>::uniqueData data(prams);

				CollisionOperator<3, 19, BGKEquilibrium, BGKCollision> BGK_BGKEQ_D3Q19;

				BGK_BGKEQ_D3Q19.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, prams,
						data);
				hasCollided = 1;
				break;
			}
			break;
		case BGK_REGULARIZED:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<3, 19> prams(stencil.getScaling(),
						viscosity, stencil, stencil.getSpeedOfSoundSquare(),
						delta_t);
				Regularized<3, 19, BGKEquilibrium>::uniqueData data(prams);
				CollisionOperator<3, 19, BGKEquilibrium, Regularized> REG_BGKEQ_D3Q19;
				REG_BGKEQ_D3Q19.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, prams,
						data);
				hasCollided = 1;
				break;
			}
			break;
		}
		break;
	case Stencil_D3Q27:
		switch (configuration->getCollisionScheme()) {
		case BGK_STANDARD:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<3, 27> prams(stencil.getScaling(),
						viscosity, stencil, stencil.getSpeedOfSoundSquare(),
						delta_t);

				CollisionOperator<3, 27, BGKEquilibrium, BGKCollision> BGK_BGKEQ_D3Q27;
				BGKCollision<3, 27, BGKEquilibrium>::uniqueData data(prams);
				BGK_BGKEQ_D3Q27.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, prams,
						data);
				hasCollided = 1;
				break;
			}
			break;

		case BGK_REGULARIZED:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				CollisionParameters<3, 27> prams(stencil.getScaling(),
						viscosity, stencil, stencil.getSpeedOfSoundSquare(),
						delta_t);
				Regularized<3, 27, BGKEquilibrium>::uniqueData data(prams);
				CollisionOperator<3, 27, BGKEquilibrium, Regularized> REG_BGKEQ_D3Q27;
				REG_BGKEQ_D3Q27.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, prams,
						data);
				hasCollided = 1;
				break;
			}
			break;
		}

		break;
	} // getStencil

	if (hasCollided == 0)
	{
		throw CollisionException ("Collision model not implemented yet -- See CollisionSelection.h");	}

} // selectCollision
} //namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_ADVANCED_COLLISIONSELECTION_H_ */
