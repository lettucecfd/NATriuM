/*
 * ErrorStats.cpp
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#include <solver/ErrorStats.h>

#include <solver/BenchmarkCFDSolver.h>

namespace natrium {
/*
template<> void ErrorStats<2>::update() {
	// this function must not be called more often than once per iteration
	// as the data for the analytic solution is constantly overwritten
	// therefor check a marker value that is set by this function (see below)

	// get analytic and numeric values
	// TODO: only assign once (see. addAnalyticSolutionToOutput)
	m_solver->m_benchmark->getAllAnalyticVelocities(m_solver->getTime(),
			m_solver->m_analyticVelocity, m_solver->m_supportPoints);
	m_solver->m_benchmark->getAllAnalyticDensities(m_solver->getTime(),
			m_solver->m_analyticDensity, m_solver->m_supportPoints);
	const vector<distributed_vector>& numericVelocity = m_solver->getVelocity();
	const distributed_vector& numericDensity = m_solver->getDensity();

	//#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
	m_maxUAnalytic = m_solver->getMaxVelocityNorm();
	// TODO take the analytic!!!
	m_solver->m_analyticDensity.add(-1.0, numericDensity);
	m_maxDensityError = m_solver->m_analyticDensity.linfty_norm();
	m_l2DensityError = m_solver->m_analyticDensity.l2_norm();
	// substract numeric from analytic velocity
	m_solver->m_analyticVelocity.at(0).add(-1.0, numericVelocity.at(0));
	m_solver->m_analyticVelocity.at(1).add(-1.0, numericVelocity.at(1));
	// calculate squares
	m_solver->m_analyticVelocity.at(0).scale(
			m_solver->m_analyticVelocity.at(0));
	m_solver->m_analyticVelocity.at(1).scale(
			m_solver->m_analyticVelocity.at(1));
	// calculate ||error (pointwise)||^2
	m_solver->m_analyticVelocity.at(0).add(m_solver->m_analyticVelocity.at(1));
	// calculate || error (pointwise) ||
	for (size_t i = 0; i < m_solver->getNumberOfDoFs(); i++) {
		m_solver->m_analyticVelocity.at(0)(i) = sqrt(
				m_solver->m_analyticVelocity.at(0)(i));
	}
	m_maxVelocityError = m_solver->m_analyticVelocity.at(0).linfty_norm();
	m_l2VelocityError = m_solver->m_analyticVelocity.at(0).l2_norm();

	// set marker value that indicates that this function has already been called
	// for the present data
	m_solver->m_analyticVelocity.at(1)(0) = 31415926;
}

template<> void ErrorStats<3>::update() {
	// this function must not be called more often than once per iteration
	// as the data for the analytic solution is constantly overwritten
	// therefor check a marker value that is set by this function (see below)

	// get analytic and numeric values
	// TODO: only assign once (see. addAnalyticSolutionToOutput)
	m_solver->m_benchmark->getAllAnalyticVelocities(m_solver->getTime(),
			m_solver->m_analyticVelocity, m_solver->m_supportPoints);
	m_solver->m_benchmark->getAllAnalyticDensities(m_solver->getTime(),
			m_solver->m_analyticDensity, m_solver->m_supportPoints);
	const vector<distributed_vector>& numericVelocity = m_solver->getVelocity();
	const distributed_vector& numericDensity = m_solver->getDensity();

	//#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
	m_maxUAnalytic = m_solver->getMaxVelocityNorm();
	// TODO take the analytic!!!
	m_solver->m_analyticDensity.add(-1.0, numericDensity);
	m_maxDensityError = m_solver->m_analyticDensity.linfty_norm();
	m_l2DensityError = m_solver->m_analyticDensity.l2_norm();
	// substract numeric from analytic velocity
	m_solver->m_analyticVelocity.at(0).add(-1.0, numericVelocity.at(0));
	m_solver->m_analyticVelocity.at(1).add(-1.0, numericVelocity.at(1));
	m_solver->m_analyticVelocity.at(2).add(-1.0, numericVelocity.at(2));
	// calculate squares
	m_solver->m_analyticVelocity.at(0).scale(
			m_solver->m_analyticVelocity.at(0));
	m_solver->m_analyticVelocity.at(1).scale(
			m_solver->m_analyticVelocity.at(1));
	m_solver->m_analyticVelocity.at(2).scale(
			m_solver->m_analyticVelocity.at(2));
	// calculate ||error (pointwise)||^2
	m_solver->m_analyticVelocity.at(0).add(m_solver->m_analyticVelocity.at(1));
	m_solver->m_analyticVelocity.at(0).add(m_solver->m_analyticVelocity.at(2));
	// calculate || error (pointwise) ||
	for (size_t i = 0; i < m_solver->getNumberOfDoFs(); i++) {
		m_solver->m_analyticVelocity.at(0)(i) = sqrt(
				m_solver->m_analyticVelocity.at(0)(i));
	}
	m_maxVelocityError = m_solver->m_analyticVelocity.at(0).linfty_norm();
	m_l2VelocityError = m_solver->m_analyticVelocity.at(0).l2_norm();

	// set marker value that indicates that this function has already been called
	// for the present data
	m_solver->m_analyticVelocity.at(1)(0) = 31415926;
} update*/

} /* namespace natrium */
