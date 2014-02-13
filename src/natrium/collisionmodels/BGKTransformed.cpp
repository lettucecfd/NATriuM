/**
 * @file BGKTransformed.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGKTransformed.h"

namespace natrium {

/// constructor
BGKTransformed::BGKTransformed(double relaxationParameter,
		boost::shared_ptr<BoltzmannModel> boltzmannModel) :
		CollisionModel(boltzmannModel), m_relaxationParameter(relaxationParameter), m_prefactor(
				-1. / (relaxationParameter + 0.5)) {


} // constructor

// destructor
BGKTransformed::~BGKTransformed() {
} // destructor

// collide
/* @note: this function is not used in the code, because it would require
 * to create the distributions vector for each dof
 */
void BGKTransformed::collideSinglePoint(vector<double>& distributions) const {

	// assert
	assert(distributions.size() == m_q);

	// calculate macroscopic entities
	double rho = m_boltzmannModel->calculateDensity(distributions);
	numeric_vector u(m_d);
	m_boltzmannModel->calculateVelocity(distributions, rho, u);

	// calculate equilibrium distribution (feq)
	vector<double> feq(m_q);
	m_boltzmannModel->getEquilibriumDistributions(feq, u, rho);

	// update distribution
	for (size_t i = 0; i < m_q; i++) {
		distributions.at(i) += m_prefactor * (distributions.at(i) - feq.at(i));
	}

} //collide

void BGKTransformed::collideSingleDoF(size_t doF, const vector<double>& feq,
		vector<distributed_vector>& f) const {
	for (size_t j = 0; j < m_q; j++) {
		f.at(j)(doF) += m_prefactor * (f.at(j)(doF) - feq.at(j));
	}

} /* collide */

}
/* namespace natrium */

