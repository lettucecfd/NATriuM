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
			const bool inInitializationProcedure, CollisionParameters<T_D, T_Q> params) const {
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());


		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;

			// Variable that stores the local distribution function values of every node
			double fLocal[T_Q];

			// Copy the needed global distribution function values into the local variable
			copyGlobalToLocalF<T_Q>(fLocal, f, i); // done

			//Calculate the local density and store it into the Parameter Handling System
			params.density = calculateDensity<T_Q>(fLocal); // done

			//Write the local density to the global density vector
			densities[i] = params.density; // write local density to global density vector

			//Calculate the local velocity and store it into the Parameter Handling System
			calculateVelocity<T_D, T_Q>(fLocal, params.velocity, params.scaling, params.density); // TODO velocities for other stencils

			//Write the local density to the global velocity matrix
			if (not inInitializationProcedure) {

				velocities.at(0)(i) = params.velocity[0] * params.scaling;
				velocities.at(1)(i) = params.velocity[1] * params.scaling;
			}

			//applyForces<T_Q>(fLocal); // TODO

			//Initialize an object of the desired collision scheme and run the relaxation process
			T_collision<T_D, T_Q, T_equilibrium> collisionScheme;
			collisionScheme.relax(fLocal, params);

			//reApplyForces<T_Q>(fLocal); // TODO

			//Finally copy the updated distribution function back to the global distribution function
			copyLocalToGlobalF<T_Q>(fLocal, f, i);



		}

	}
};



} /* namespace natrium */

#endif /* COLLISIONOPERATOR_H_ */
