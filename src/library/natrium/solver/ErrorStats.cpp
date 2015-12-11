/*
 * ErrorStats.cpp
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#include "ErrorStats.h"

#include "BenchmarkCFDSolver.h"

namespace natrium {

template<size_t dim>
ErrorStats<dim>::ErrorStats(BenchmarkCFDSolver<dim> * cfdsolver,
		const std::string tableFileName) :
		m_solver(cfdsolver), m_filename(tableFileName), m_outputOff(
				tableFileName == "") {
	// set information
	m_iterationNumber = 100000000037;
	m_time = 0.0;
	m_maxVelocityError = 0.0;
	m_maxDensityError = 0.0;
	m_l2VelocityError = 0.0;
	m_l2DensityError = 0.0;
	m_maxUAnalytic = 0.0;

	// create file (if necessary)
	if (m_solver->getIterationStart() > 0) {
		m_errorsTableFile = boost::make_shared<std::fstream>(tableFileName,
				std::fstream::out | std::fstream::app);
	} else {
		m_errorsTableFile = boost::make_shared<std::fstream>(tableFileName,
				std::fstream::out);
		printHeaderLine();
	}
	MPI_sync();
}
template ErrorStats<2>::ErrorStats(BenchmarkCFDSolver<2> * cfdsolver,
		const std::string tableFileName);
template ErrorStats<3>::ErrorStats(BenchmarkCFDSolver<3> * cfdsolver,
		const std::string tableFileName);

template<size_t dim>
void ErrorStats<dim>::update() {
	// this function must not be called more often than once per iteration
	// as the data for the analytic solution is constantly overwritten
	// therefor check a marker value that is set by this function (see below)
	if (m_iterationNumber == m_solver->getIteration()) {
		return;
	}
	m_iterationNumber = m_solver->getIteration();
	m_time = m_solver->getTime();
	// get analytic and numeric values
	// TODO: only assign once (see. addAnalyticSolutionToOutput)
	m_solver->getAllAnalyticVelocities(m_solver->getTime(),
			m_solver->m_analyticVelocity, m_solver->m_supportPoints);
	m_solver->getAllAnalyticDensities(m_solver->getTime(),
			m_solver->m_analyticDensity, m_solver->m_supportPoints);
	const vector<distributed_vector>& numericVelocity = m_solver->getVelocity();
	const distributed_vector& numericDensity = m_solver->getDensity();

	//#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
	m_solver->m_analyticDensity.add(-1.0, numericDensity);
	m_maxDensityError = m_solver->m_analyticDensity.linfty_norm();
	m_l2DensityError = m_solver->m_analyticDensity.l2_norm();

	// calculate maximum analytic velocity norm
	const dealii::IndexSet& locally_owned_dofs = m_solver->getAdvectionOperator()->getLocallyOwnedDofs();
	m_maxUAnalytic = Math::maxVelocityNorm(m_solver->m_analyticVelocity, locally_owned_dofs);
	m_l2UAnalytic = Math::velocity2Norm(m_solver->m_analyticVelocity, locally_owned_dofs);

	// substract numeric from analytic velocity
	m_solver->m_analyticVelocity.at(0).add(-1.0, numericVelocity.at(0));
	m_solver->m_analyticVelocity.at(1).add(-1.0, numericVelocity.at(1));
	if (dim == 3) {
		m_solver->m_analyticVelocity.at(2).add(-1.0, numericVelocity.at(2));
	}
	// calculate squares
	m_solver->m_analyticVelocity.at(0).scale(
			m_solver->m_analyticVelocity.at(0));
	m_solver->m_analyticVelocity.at(1).scale(
			m_solver->m_analyticVelocity.at(1));
	if (dim == 3) {
		m_solver->m_analyticVelocity.at(2).scale(
				m_solver->m_analyticVelocity.at(2));
	}
	// calculate ||error (pointwise)||^2
	m_solver->m_analyticVelocity.at(0).add(m_solver->m_analyticVelocity.at(1));
	if (dim == 3) {
		m_solver->m_analyticVelocity.at(0).add(
				m_solver->m_analyticVelocity.at(2));
	}
	// calculate || error (pointwise) ||
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++){
		size_t i = *it;
		m_solver->m_analyticVelocity.at(0)(i) = sqrt(
				m_solver->m_analyticVelocity.at(0)(i));
	}
	m_maxVelocityError = m_solver->m_analyticVelocity.at(0).linfty_norm();
	m_l2VelocityError = m_solver->m_analyticVelocity.at(0).l2_norm();

	// set marker value that indicates that this function has already been called
	// for the present data
	m_solver->m_analyticVelocity.at(1)(0) = 31415926;
} /*update*/
template void ErrorStats<2>::update();
template void ErrorStats<3>::update();

template<size_t dim>
void ErrorStats<dim>::printHeaderLine() {
	if (is_MPI_rank_0()) {
		(*m_errorsTableFile)
				<< "#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
				<< endl;
	}
}
template void ErrorStats<2>::printHeaderLine();
template void ErrorStats<3>::printHeaderLine();

template<size_t dim>
void ErrorStats<dim>::printNewLine() {
	if (not isUpToDate()) {
		update();
	}
	if (is_MPI_rank_0()) {
		(*m_errorsTableFile) << m_iterationNumber << " " << m_time << " "
				<< m_maxUAnalytic << " " << m_maxVelocityError << " "
				<< m_maxDensityError << " " << m_l2VelocityError << " "
				<< m_l2DensityError << endl;
	}
}
template void ErrorStats<2>::printNewLine();
template void ErrorStats<3>::printNewLine();

} /* namespace natrium */
