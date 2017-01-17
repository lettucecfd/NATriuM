/*
 * CollisionOperator.h
 *
 *  Created on: 06.01.2017
 *      Author: dominik
 */

#ifndef COLLISIONOPERATOR_H_
#define COLLISIONOPERATOR_H_

#include <vector>

namespace natrium {

class CollisionOperator {
public:
	CollisionOperator(boost::shared_ptr<SolverConfiguration>& configuration);
	virtual ~CollisionOperator();

	void collide(boost::shared_ptr<SolverConfiguration>& configuration, DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	template<int T_D, int T_Q , class T_collision, class T_equilibrium>
	void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;


};

} /* namespace natrium */

#endif /* COLLISIONOPERATOR_H_ */
