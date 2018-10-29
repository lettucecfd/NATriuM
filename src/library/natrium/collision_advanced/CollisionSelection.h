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

/**
 * This macro is only used in collision selection, to avoid the repetition for each equilibrium, collision, etc.
 * The macro performs the collision step.
 * It is important at this point that everything is specialized at runtime so that the compiler is able to inline
 * the collision functions.
 * This circumvents the use of inheritance, memory allocation and other expensive operations in the collision loop.
 */
#define DEFINE_POSSIBLE_COLLISION(T_D, T_Q, STENCIL, COLLISION, EQUILIBRIUM, collision_type, equilibrium_type) \
		if ( (configuration.getStencil() == STENCIL) && \
				(configuration.getCollisionScheme() == COLLISION) && \
				(configuration.getEquilibriumScheme() == EQUILIBRIUM) ){ \
				GeneralCollisionData<T_D, T_Q> genData(configuration, problemDescription, \
						stencil.getScaling(), viscosity, stencil, \
						stencil.getSpeedOfSoundSquare(), delta_t); \
				typename collision_type<T_D, T_Q, equilibrium_type>::SpecificCollisionData\
					specData(genData);\
				CollisionOperator<T_D, T_Q, equilibrium_type, collision_type> collision_operator;\
				collision_operator.collideAll(f, densities, velocities, locally_owned_dofs,\
						inInitializationProcedure, genData, specData); \
				has_collided += 1; \
		}

#define DEFINE_POSSIBLE_COMPRESSIBLE_COLLISION(T_D, T_Q, STENCIL, COLLISION, EQUILIBRIUM, collision_type, equilibrium_type) \
		if ( (configuration.getStencil() == STENCIL) && \
				(configuration.getCollisionScheme() == COLLISION) && \
				(configuration.getEquilibriumScheme() == EQUILIBRIUM) ){ \
				GeneralCollisionData<T_D, T_Q> genData(configuration, problemDescription, \
						stencil.getScaling(), viscosity, stencil, \
						stencil.getSpeedOfSoundSquare(), delta_t); \
				typename collision_type<T_D, T_Q, equilibrium_type>::SpecificCollisionData\
					specData(genData);\
				CollisionOperator<T_D, T_Q, equilibrium_type, collision_type> collision_operator;\
				collision_operator.collideAll(f, densities, velocities, locally_owned_dofs,\
						inInitializationProcedure, genData, specData); \
				has_collided += 1; \
		}


/**
 *  @short select and perform collision
 */
template<size_t T_D>
inline void selectCollision(SolverConfiguration& configuration,
		const ProblemDescription<T_D>& problemDescription,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure);

template<>
inline void selectCollision<2>(SolverConfiguration& configuration,
		const ProblemDescription<2>& problemDescription,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure) {
	// CAUTION: The names of the arguments are used in the macro

	// CAUTION: The name has_collided is used in the macro
	int has_collided = false;

	// Select and perform collision (the macro is defined above)
	// Caution: You need to make sure that no configuration appears twice here.
	// This would cause multiple collision steps per time step
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, BGK_REGULARIZED, BGK_EQUILIBRIUM, Regularized, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, MRT_STANDARD, BGK_EQUILIBRIUM, MultipleRelaxationTime, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,25, Stencil_D2Q25, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,25, Stencil_D2Q25H, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, BGK_STANDARD, INCOMPRESSIBLE_EQUILIBRIUM, BGKCollision, IncompressibleEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, BGK_REGULARIZED, INCOMPRESSIBLE_EQUILIBRIUM, Regularized, IncompressibleEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, MRT_STANDARD, INCOMPRESSIBLE_EQUILIBRIUM, MultipleRelaxationTime, IncompressibleEquilibrium);

	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, BGK_STANDARD, STEADYSTATE_EQUILIBRIUM, BGKCollision, SteadyStateEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, BGK_REGULARIZED, STEADYSTATE_EQUILIBRIUM, Regularized, SteadyStateEquilibrium);
	DEFINE_POSSIBLE_COLLISION(2,9, Stencil_D2Q9, MRT_STANDARD, STEADYSTATE_EQUILIBRIUM, MultipleRelaxationTime, SteadyStateEquilibrium);

	if (not has_collided) {
		throw CollisionException(
				"Severe error: Collision model not implemented yet -- cf. CollisionSelection.h");
	}
	if (has_collided > 1) {
		throw CollisionException(
				"Severe error: More than one collision model was executed per time step. "
				"This means that there is a bug in CollisionSelection.h");

	}

} // selectCollision 2D

template<size_t T_D>
inline void selectCollision(SolverConfiguration& configuration,
		const ProblemDescription<T_D>& problemDescription,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities, distributed_vector& temperature,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure);

template<>
inline void selectCollision<2>(SolverConfiguration& configuration,
		const ProblemDescription<2>& problemDescription,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities, distributed_vector& temperature,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure) {
	// CAUTION: The names of the arguments are used in the macro

	// CAUTION: The name has_collided is used in the macro
	int has_collided = false;

	// Select and perform collision (the macro is defined above)
	// Caution: You need to make sure that no configuration appears twice here.
	// This would cause multiple collision steps per time step

	DEFINE_POSSIBLE_COMPRESSIBLE_COLLISION(2,25, Stencil_D2Q25H, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);

	if (not has_collided) {
		throw CollisionException(
				"Severe error: Collision model not implemented yet -- cf. CollisionSelection.h");
	}
	if (has_collided > 1) {
		throw CollisionException(
				"Severe error: More than one collision model was executed per time step. "
				"This means that there is a bug in CollisionSelection.h");

	}

} // selectCollision 2D



template<>
inline void selectCollision<3>(SolverConfiguration& configuration,
		const ProblemDescription<3>& problemDescription,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs, const double viscosity,
		const double delta_t, const Stencil& stencil,
		const bool inInitializationProcedure) {
	// CAUTION: The names of the arguments are used in the macro

	// CAUTION: The name has_collided is used in the macro
	bool has_collided = false;

	DEFINE_POSSIBLE_COLLISION(3,15, Stencil_D3Q15, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,15, Stencil_D3Q15, BGK_REGULARIZED, BGK_EQUILIBRIUM, Regularized, BGKEquilibrium);

	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, BGK_REGULARIZED, BGK_EQUILIBRIUM, Regularized, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, MRT_STANDARD, BGK_EQUILIBRIUM, MultipleRelaxationTime, BGKEquilibrium);

	DEFINE_POSSIBLE_COLLISION(3,27, Stencil_D3Q27, BGK_STANDARD, BGK_EQUILIBRIUM, BGKCollision, BGKEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,27, Stencil_D3Q27, BGK_REGULARIZED, BGK_EQUILIBRIUM, Regularized, BGKEquilibrium);
	// ---------
	DEFINE_POSSIBLE_COLLISION(3,15, Stencil_D3Q15, BGK_STANDARD, STEADYSTATE_EQUILIBRIUM, BGKCollision, SteadyStateEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,15, Stencil_D3Q15, BGK_REGULARIZED, STEADYSTATE_EQUILIBRIUM, Regularized, SteadyStateEquilibrium);

	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, BGK_STANDARD, STEADYSTATE_EQUILIBRIUM, BGKCollision, SteadyStateEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, BGK_REGULARIZED, STEADYSTATE_EQUILIBRIUM, Regularized, SteadyStateEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, MRT_STANDARD, STEADYSTATE_EQUILIBRIUM, MultipleRelaxationTime, SteadyStateEquilibrium);

	DEFINE_POSSIBLE_COLLISION(3,27, Stencil_D3Q27, BGK_STANDARD, STEADYSTATE_EQUILIBRIUM, BGKCollision, SteadyStateEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,27, Stencil_D3Q27, BGK_REGULARIZED, STEADYSTATE_EQUILIBRIUM, Regularized, SteadyStateEquilibrium);
	// -------------
	DEFINE_POSSIBLE_COLLISION(3,15, Stencil_D3Q15, BGK_STANDARD, INCOMPRESSIBLE_EQUILIBRIUM, BGKCollision, IncompressibleEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,15, Stencil_D3Q15, BGK_REGULARIZED, INCOMPRESSIBLE_EQUILIBRIUM, Regularized, IncompressibleEquilibrium);

	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, BGK_STANDARD, INCOMPRESSIBLE_EQUILIBRIUM, BGKCollision, IncompressibleEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, BGK_REGULARIZED, INCOMPRESSIBLE_EQUILIBRIUM, Regularized, IncompressibleEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,19, Stencil_D3Q19, MRT_STANDARD, INCOMPRESSIBLE_EQUILIBRIUM, MultipleRelaxationTime, IncompressibleEquilibrium);

	DEFINE_POSSIBLE_COLLISION(3,27, Stencil_D3Q27, BGK_STANDARD, INCOMPRESSIBLE_EQUILIBRIUM, BGKCollision, IncompressibleEquilibrium);
	DEFINE_POSSIBLE_COLLISION(3,27, Stencil_D3Q27, BGK_REGULARIZED, INCOMPRESSIBLE_EQUILIBRIUM, Regularized, IncompressibleEquilibrium);

	if (not has_collided) {
		throw CollisionException(
				"Collision model not implemented yet -- cf. CollisionSelection.h");
	}

} // selectCollision 3D
} //namespace natrium

#endif /* LIBRARY_NATRIUM_COLLISION_ADVANCED_COLLISIONSELECTION_H_ */
