/*
 * CollisionOperator.h
 *
 *  Created on: 06.01.2017
 *      Author: dominik
 */

#ifndef COLLISIONOPERATOR_H_
#define COLLISIONOPERATOR_H_

#include "../stencils/Stencil.h"
#include "AuxiliaryCollisionFunctions.h"
#include <vector>
#include "Equilibria.h"
#include "CollisionSchemes.h"



namespace natrium {




template<int T_D, int T_Q, template <int T_D, int T_Q> class T_equilibrium, template <int T_D, int T_Q, template <int T_D, int T_Q> class T_equilibrium> class T_collision>
class CollisionOperator {
public:
/*
	void collide(boost::shared_ptr<SolverConfiguration>& configuration, DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const; */


	void collideAll(DistributionFunctions& f,
			distributed_vector& densities, vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const {
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;

			CollisionParameters<T_D,T_Q> params;
			params.scaling = 1.0; // TODO adjust to global scaling
			params.tau = 1.0 // TODO
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

			T_collision<T_D, T_Q, T_equilibrium> collisionScheme;
			collisionScheme.relax(fLocal, params);  //TODO

			copyLocalToGlobalF<T_Q>(fLocal,f,i);

			//reApplyForces<T_Q>(fLocal); // TODO


		}

	}
};


inline void selectCollide(
		boost::shared_ptr<SolverConfiguration>& configuration,
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) {

	switch (configuration->getStencil()) {
	case Stencil_D2Q9:
		switch (configuration->getCollisionScheme()) {
		case BGK_STANDARD:
			switch (configuration->getEquilibriumScheme()) {
			case BGK_EQUILIBRIUM:
				cout << "ok";

			CollisionOperator<2,9,BGKEquilibrium, BGKCollision> test;
			test.collideAll (f, densities,velocities, locally_owned_dofs, inInitializationProcedure);
			}
		}

	}
}


} /* namespace natrium */

#endif /* COLLISIONOPERATOR_H_ */
