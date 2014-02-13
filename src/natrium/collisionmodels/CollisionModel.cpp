/**
 * @file CollisionModel.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CollisionModel.h"

#include <cmath>

namespace natrium {

// constructor
CollisionModel::CollisionModel(boost::shared_ptr<BoltzmannModel> boltzmannModel) :
		m_boltzmannModel(boltzmannModel), m_d(boltzmannModel->getD()), m_q(
				boltzmannModel->getQ()) {

} // constructor

CollisionModel::~CollisionModel() {
}

void CollisionModel::collideAll(vector<distributed_vector>& f,
		distributed_vector& densities,
		vector<distributed_vector>& velocities) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = m_boltzmannModel->getQ();

	assert(f.size() == Q);
	assert(velocities.size() == m_d);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

//for all degrees of freedom
	for (size_t i = 0; i < n_dofs; i++) {

		// calculate density
		densities(i) = 0;
		for (size_t j = 0; j < Q; j++) {
			densities(i) += f.at(j)(i);
		}

		// calculate velocity
		// for all velocity components
		for (size_t j = 0; j < m_d; j++) {
			velocities.at(j)(i) = 0;
			for (size_t k = 0; k < Q; k++) {
				// check that the components of the direction vectors are in fact ints
#ifdef DEBUG
				double integer;
				assert(
						std::modf(m_boltzmannModel->getDirection(k)(j), &integer)
						== 0.0);
#endif
				switch (int(m_boltzmannModel->getDirection(k)(j))) {
				case 0:
					break;
				case 1:
					velocities.at(j)(i) += f.at(k)(i);
					break;
				case -1:
					velocities.at(j)(i) -= f.at(k)(i);
					break;
				default:
					velocities.at(j)(i) += f.at(k)(i)
							* m_boltzmannModel->getDirection(k)(j);
				}
			}
			assert(densities(i) > 1e-10);
			velocities.at(j)(i) /= densities(i);
		}

		// calculate equilibrium distribution
		// TODO Optimize by passing different arguments
		vector<double> feq(Q, 0.0);
		numeric_vector u(m_d);
		for (size_t j = 0; j < m_d; j++) {
			u(j) = velocities.at(j)(i);
		}
		m_boltzmannModel->getEquilibriumDistributions(feq, u, densities(i));

		// BGK collision
		collideSingleDoF(i, feq, f);
	}
} /*collideAll */

} /* namespace natrium */
