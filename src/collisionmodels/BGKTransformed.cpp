/**
 * @file BGKTransformed.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGKTransformed.h"

namespace natrium {

/// constructor
BGKTransformed::BGKTransformed(float_t relaxationParameter,
		boost::shared_ptr<BoltzmannModel> boltzmannModel):
		CollisionModel(relaxationParameter, boltzmannModel),
		m_prefactor(-1./(relaxationParameter + 0.5)){

} // constructor

// destructor
BGKTransformed::~BGKTransformed() {
} // destructor


// collide
void BGKTransformed::collide(vector<float_t>& distributions) const {

	// assert
	assert(distributions.size() == m_q);

	// calculate macroscopic entities
	float_t rho = m_boltzmannModel->calculateDensity(distributions);
	numeric_vector u(m_d);
	m_boltzmannModel->calculateVelocity(distributions, rho, u);

	// calculate equilibrium distribution (feq)
	vector<float_t> feq(m_q);
	m_boltzmannModel->getEquilibriumDistributions(feq, u, rho);

	// update distribution
	for (size_t i = 0; i < m_q; i++){
		distributions.at(i) += m_prefactor * (distributions.at(i) - feq.at(i));
	}

} //collide


} /* namespace natrium */
