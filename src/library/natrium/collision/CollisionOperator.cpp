/*
 * CollisionOperator.cpp
 *
 *  Created on: 06.01.2017
 *      Author: dominik
 */

#include "CollisionOperator.h"
#include "../stencils/Stencil.h"
#include "AuxiliaryCollisionFunctions.h"

namespace natrium {

CollisionOperator::CollisionOperator(boost::shared_ptr<SolverConfiguration>& configuration) {
	// TODO Auto-generated constructor stub

}

CollisionOperator::~CollisionOperator() {
	// TODO Auto-generated destructor stub
}

 inline void CollisionOperator::collide(boost::shared_ptr<SolverConfiguration>& configuration, DistributionFunctions& f,
		distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {



switch(configuration->getStencil()){
case Stencil_D2Q9:
	switch(configuration->getCollisionScheme())
	{
	case BGK_STANDARD:
		collideAll<2,9,BGKStandard,BGKStandard>(f, densities, velocities, locally_owned_dofs, inInitializationProcedure);
	}

}
 }

template<int T_D, int T_Q, class T_collision, class T_equilibrium>
void CollisionOperator::collideAll(DistributionFunctions& f,
		distributed_vector& densities,
		vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const
{
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		double density;
		double velocity[T_Q];
		double fLocal[T_Q];
		double feqLocal[T_Q];

		copyGlobalToLocalF<T_Q>(fLocal,f,i); // done

		calculateDensity<T_Q>(fLocal,density); // done

		calculateVelocity<T_D,T_Q>(fLocal,velocity); // TODO

		applyForces<T_Q>(fLocal); // TODO

		relaxDistributions<T_D,T_Q,T_collision,T_equilibrium>(fLocal,feqLocal);  //TODO

		reApplyForces<T_Q>(fLocal); // TODO






}
}
}/* namespace natrium */
