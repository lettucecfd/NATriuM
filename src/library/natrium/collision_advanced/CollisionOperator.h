/*
 * CollisionOperator.h
 *
 *  Created on: 06.01.2017
 *      Author: dominik
 */

#ifndef COLLISIONOPERATOR_H_
#define COLLISIONOPERATOR_H_

#include <vector>
#include "../stencils/Stencil.h"
#include "AuxiliaryCollisionFunctions.h"
#include "Equilibria.h"
#include "CollisionSchemes.h"

namespace natrium {

template<int T_D, int T_Q, template<int T_D, int T_Q> class T_equilibrium,
		template<int T_D, int T_Q,
				template<int T_D, int T_Q> class T_equilibrium> class T_collision>
class CollisionOperator {
public:

	void collideAll(DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			const bool inInitializationProcedure,
			GeneralCollisionData<T_D, T_Q>& genData,
			typename T_collision<T_D, T_Q, T_equilibrium>::SpecificCollisionData& specData) const {
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

		T_collision<T_D, T_Q, T_equilibrium> collisionScheme;

		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;

			// Variable that stores the local distribution function values of every node

			// Copy the needed global distribution function values into the local variable
			copyGlobalToLocalF<T_Q>(genData.fLocal, f, i); // done

			//Calculate the local density and store it into the Parameter Handling System
			genData.density = calculateDensity<T_Q>(genData.fLocal); // done
			if (1e-10 >= genData.density) {
				throw DensityZeroException(
						"Density too small in collision. Decrease time step size.");
			}

			//Write the local density to the global density vector
			densities[i] = genData.density; // write local density to global density vector

			//Calculate the local velocity and store it into the Parameter Handling System
			calculateVelocity<T_D, T_Q>(genData.fLocal, genData.velocity,
					genData.density, genData); // TODO velocities for other stencils

			//Write the local density to the global velocity matrix
			if (not inInitializationProcedure) {
				for (size_t j = 0; j < T_D; ++j) {
					velocities.at(j)(i) = genData.velocity[j] * genData.scaling;
				}
				if (genData.problemDescription.hasExternalForce()) {
					applyMacroscopicForces<T_D, T_Q>(velocities, i, genData);
					applyForces<T_D, T_Q>(genData);
				}
			} else {
				for (size_t j = 0; j < T_D; ++j) {
					genData.velocity[i] = velocities.at(0)(i);
				}
			}

			//Initialize an object of the desired collision scheme and run the relaxation process

			collisionScheme.relax(genData.fLocal, genData, specData);

			//reApplyForces<T_Q>(fLocal); // TODO

			//Finally copy the updated distribution function back to the global distribution function
			copyLocalToGlobalF<T_Q>(genData.fLocal, f, i);

		}

	}
};

} /* namespace natrium */

#endif /* COLLISIONOPERATOR_H_ */