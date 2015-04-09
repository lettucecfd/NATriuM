
/**
 * @file CollisionModel.cpp
 * @short
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CollisionModel.h"

namespace natrium {

void CollisionModel::getEquilibriumDistributions(vector<double>& feq,
			const numeric_vector& u, const double rho) const{
	assert (feq.size() == getStencil()->getQ());
	for (size_t i = 0; i < getStencil()->getQ(); i++){
		feq.at(i) = getEquilibriumDistribution( i, u, rho);
	}
}

}
