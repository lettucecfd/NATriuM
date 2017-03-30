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
			bool inInitializationProcedure, CollisionParameters<T_D, T_Q> params) const {
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

	//	CollisionParameters<T_D, T_Q> params;

		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;


			//params.scaling = 1.0; // TODO adjust to global scaling
			//params.tau = 3*params.viscosity+0.5; // TODO
			double fLocal[T_Q];

			copyGlobalToLocalF<T_Q>(fLocal, f, i); // done

			params.density = calculateDensity<T_Q>(fLocal); // done

			densities[i] = params.density; // write local density to global density vector

			calculateVelocity<T_D, T_Q>(fLocal, params.velocity, 1.0, 1.0); // TODO velocities for other stencils

			if (not inInitializationProcedure) {

				velocities.at(0)(i) = params.velocity[0] * params.scaling;
				velocities.at(1)(i) = params.velocity[1] * params.scaling;
			}

			//applyForces<T_Q>(fLocal); // TODO

			T_collision<T_D, T_Q, T_equilibrium> collisionScheme(params);
			collisionScheme.relax(fLocal, params);  //TODO

			//reApplyForces<T_Q>(fLocal); // TODO

			copyLocalToGlobalF<T_Q>(fLocal, f, i); //TODO



		}

	}
};



} /* namespace natrium */

#endif /* COLLISIONOPERATOR_H_ */
