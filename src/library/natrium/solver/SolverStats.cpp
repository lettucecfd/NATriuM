/*
 * SolverStats.cpp
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#include <solver/SolverStats.h>

namespace natrium {

/*

template<>
void natrium::SolverStats<2>::calulateResiduals(size_t iteration) {

	assert(iteration % 10 == 0);
	if (not (m_solver->m_i - m_solver->m_iterationStart < 100)) {
		// i.e. not first visit of this if-statement
		///// CALCULATE MAX VELOCITY VARIATION
		// substract new from old velocity
		m_solver->m_tmpVelocity.at(0).add(-1.0, m_solver->m_velocity.at(0));
		m_solver->m_tmpVelocity.at(1).add(-1.0, m_solver->m_velocity.at(1));
		// calculate squares
		m_solver->m_tmpVelocity.at(0).scale(m_solver->m_tmpVelocity.at(0));
		m_solver->m_tmpVelocity.at(1).scale(m_solver->m_tmpVelocity.at(1));
		// calculate ||error (pointwise)||^2
		m_solver->m_tmpVelocity.at(0).add(m_solver->m_tmpVelocity.at(1));
		m_solver->m_residuumVelocity = sqrt(
				m_solver->m_tmpVelocity.at(0).linfty_norm());
		///// CALCULATE MAX DENSITY VARIATION
		m_solver->m_tmpDensity.add(-1.0, m_solver->m_density);
		m_solver->m_residuumDensity = sqrt(
				m_solver->m_tmpVelocity.at(0).linfty_norm());

	}
	// (if first visit of this if-statement, only this part is executed)
	assert(m_solver->m_tmpVelocity.size() == m_solver->m_velocity.size());
	m_solver->m_tmpVelocity.at(0) = m_solver->m_velocity.at(0);
	m_solver->m_tmpVelocity.at(1) = m_solver->m_velocity.at(1);
	m_solver->m_tmpDensity = m_solver->m_density;
}
template<>
void natrium::SolverStats<3>::calulateResiduals(size_t iteration) {

	assert(iteration % 10 == 0);
	if (not (m_solver->m_i - m_solver->m_iterationStart < 100)) {
		// i.e. not first visit of this if-statement
		///// CALCULATE MAX VELOCITY VARIATION
		// substract new from old velocity
		m_solver->m_tmpVelocity.at(0).add(-1.0, m_solver->m_velocity.at(0));
		m_solver->m_tmpVelocity.at(1).add(-1.0, m_solver->m_velocity.at(1));
		m_solver->m_tmpVelocity.at(2).add(-1.0, m_solver->m_velocity.at(2));
		// calculate squares
		m_solver->m_tmpVelocity.at(0).scale(m_solver->m_tmpVelocity.at(0));
		m_solver->m_tmpVelocity.at(1).scale(m_solver->m_tmpVelocity.at(1));
		m_solver->m_tmpVelocity.at(2).scale(m_solver->m_tmpVelocity.at(2));
		// calculate ||error (pointwise)||^2
		m_solver->m_tmpVelocity.at(0).add(m_solver->m_tmpVelocity.at(1));
		m_solver->m_tmpVelocity.at(0).add(m_solver->m_tmpVelocity.at(2));
		m_solver->m_residuumVelocity = sqrt(
				m_solver->m_tmpVelocity.at(0).linfty_norm());
		///// CALCULATE MAX DENSITY VARIATION
		m_solver->m_tmpDensity.add(-1.0, m_solver->m_density);
		m_solver->m_residuumDensity = sqrt(
				m_solver->m_tmpVelocity.at(0).linfty_norm());

	}
	// (if first visit of this if-statement, only this part is executed)
	assert(m_solver->m_tmpVelocity.size() == m_solver->m_velocity.size());
	m_solver->m_tmpVelocity.at(0) = m_solver->m_velocity.at(0);
	m_solver->m_tmpVelocity.at(1) = m_solver->m_velocity.at(1);
	m_solver->m_tmpVelocity.at(2) = m_solver->m_velocity.at(2);

	m_solver->m_tmpDensity = m_solver->m_density;
}
*/
} /* namespace natrium */
