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

template<size_t T_D>
inline void selectCollision(const boost::shared_ptr<SolverConfiguration>& configuration,
		const boost::shared_ptr<ProblemDescription<T_D>>& problemDescription,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure);

template<>
inline void selectCollision<2>(
		const boost::shared_ptr<SolverConfiguration>& configuration,
		const boost::shared_ptr<ProblemDescription<2>>& problemDescription,
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
				GeneralCollisionData<2, 9> genData(configuration,
						problemDescription, stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);

				BGKCollision<2, 9, BGKEquilibrium>::SpecificCollisionData specData(
						genData);

				CollisionOperator<2, 9, BGKEquilibrium, BGKCollision> BGK_BGKEQ_D2Q9;

				BGK_BGKEQ_D2Q9.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, genData,
						specData);
				hasCollided = 1;
				break;
			}
			break;
		case BGK_REGULARIZED:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				GeneralCollisionData<2, 9> genData(configuration,
						problemDescription, stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);

				Regularized<2, 9, BGKEquilibrium>::SpecificCollisionData specData(
						genData);

				CollisionOperator<2, 9, BGKEquilibrium, Regularized> REG_BGKEQ_D2Q9;

				REG_BGKEQ_D2Q9.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, genData,
						specData);
				hasCollided = 1;
				break;
			}
			break;
		}
		break;

}

	if (hasCollided == 0) {
			throw CollisionException(
					"Collision model not implemented yet -- See CollisionSelection.h");
	}

} // selectCollision 2D

	template<>
	inline void selectCollision<3>(
			const boost::shared_ptr<SolverConfiguration>& configuration,
			const boost::shared_ptr<ProblemDescription<3>>& problemDescription,
			DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs, const double viscosity,
			const double delta_t, const Stencil& stencil,
			const bool inInitializationProcedure) {

		bool hasCollided = 0;

	switch(configuration->getStencil())
	{

	case Stencil_D3Q19:
		switch (configuration->getCollisionScheme()) {
		case BGK_STANDARD:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				GeneralCollisionData<3, 19> genData(configuration,
						problemDescription, stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);

				BGKCollision<3, 19, BGKEquilibrium>::SpecificCollisionData specData(
						genData);

				CollisionOperator<3, 19, BGKEquilibrium, BGKCollision> BGK_BGKEQ_D3Q19;

				BGK_BGKEQ_D3Q19.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, genData,
						specData);
				hasCollided = 1;
				break;
			}
			break;
		case BGK_REGULARIZED:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				GeneralCollisionData<3, 19> genData(configuration,
						problemDescription, stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);
				Regularized<3, 19, BGKEquilibrium>::SpecificCollisionData specData(
						genData);
				CollisionOperator<3, 19, BGKEquilibrium, Regularized> REG_BGKEQ_D3Q19;
				REG_BGKEQ_D3Q19.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, genData,
						specData);
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
				GeneralCollisionData<3, 27> genData(configuration,
						problemDescription, stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);

				CollisionOperator<3, 27, BGKEquilibrium, BGKCollision> BGK_BGKEQ_D3Q27;
				BGKCollision<3, 27, BGKEquilibrium>::SpecificCollisionData specData(
						genData);
				BGK_BGKEQ_D3Q27.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, genData,
						specData);
				hasCollided = 1;
				break;
			}
			break;

		case BGK_REGULARIZED:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				GeneralCollisionData<3, 27> genData(configuration,
						problemDescription, stencil.getScaling(), viscosity,
						stencil, stencil.getSpeedOfSoundSquare(), delta_t);
				Regularized<3, 27, BGKEquilibrium>::SpecificCollisionData specData(
						genData);
				CollisionOperator<3, 27, BGKEquilibrium, Regularized> REG_BGKEQ_D3Q27;
				REG_BGKEQ_D3Q27.collideAll(f, densities, velocities,
						locally_owned_dofs, inInitializationProcedure, genData,
						specData);
				hasCollided = 1;
				break;
			}
			break;
		}

		break;
	} // getStencil

	if (hasCollided == 0) {
		throw CollisionException(
				"Collision model not implemented yet -- See CollisionSelection.h");
	}

} // selectCollision 3D
} //namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_ADVANCED_COLLISIONSELECTION_H_ */
