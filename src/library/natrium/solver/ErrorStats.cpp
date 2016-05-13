/*
 * ErrorStats.cpp
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#include "ErrorStats.h"

#include "BenchmarkCFDSolver.h"
#include "deal.II/numerics/vector_tools.h"
#include "deal.II/base/mpi.h"
#include "math.h"

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

	// calculate maximum analytic velocity norm
	const dealii::IndexSet& locally_owned_dofs =
			m_solver->getAdvectionOperator()->getLocallyOwnedDofs();
	m_maxUAnalytic = Math::maxVelocityNorm(m_solver->m_analyticVelocity,
			locally_owned_dofs);
	m_l2UAnalytic = Math::velocity2Norm(m_solver->m_analyticVelocity,
			locally_owned_dofs);

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
	for (; it != end; it++) {
		size_t i = *it;
		m_solver->m_analyticVelocity.at(0)(i) = sqrt(
				m_solver->m_analyticVelocity.at(0)(i));
	}

	//rho
	const dealii::Function<dim> & f_rho =
			*m_solver->m_benchmark->getAnalyticRhoFunction(m_solver->getTime());
	dealii::Vector<double> local_errors(
			m_solver->getProblemDescription()->getMesh()->n_active_cells());
	dealii::VectorTools::integrate_difference(
			m_solver->getAdvectionOperator()->getMapping(),
			*m_solver->getAdvectionOperator()->getDoFHandler(), numericDensity,
			f_rho, local_errors,
			*m_solver->getAdvectionOperator()->getQuadrature(),
			dealii::VectorTools::L2_norm);
	const double local_rho_error = local_errors.l2_norm();
	const double lsq = local_rho_error * local_rho_error;
	m_l2DensityError = sqrt(dealii::Utilities::MPI::sum(lsq,
	MPI_COMM_WORLD));

	// u
	const dealii::Function<dim>& f_u =
			*m_solver->m_benchmark->getAnalyticUFunction(m_solver->getTime());
	AnalyticU<dim> ana_ux(f_u, 0);
	AnalyticU<dim> ana_uy(f_u, 1);
	// ux
	local_errors = 0;
	dealii::VectorTools::integrate_difference(
			m_solver->getAdvectionOperator()->getMapping(),
			*m_solver->getAdvectionOperator()->getDoFHandler(),
			numericVelocity.at(0), ana_ux, local_errors,
			*m_solver->getAdvectionOperator()->getQuadrature(),
			dealii::VectorTools::L2_norm);
	double one_component_local_error = local_errors.l2_norm();
	double total_local_error = one_component_local_error
			* one_component_local_error;
	// uy
	local_errors = 0;
	dealii::VectorTools::integrate_difference(
			m_solver->getAdvectionOperator()->getMapping(),
			*m_solver->getAdvectionOperator()->getDoFHandler(),
			numericVelocity.at(1), ana_uy, local_errors,
			*m_solver->getAdvectionOperator()->getQuadrature(),
			dealii::VectorTools::L2_norm);
	one_component_local_error = local_errors.l2_norm();
	total_local_error += one_component_local_error * one_component_local_error;
	// uz
	if (dim == 3) {
		AnalyticU<dim> ana_uz(f_u, 2);
		local_errors = 0;
		dealii::VectorTools::integrate_difference(
				m_solver->getAdvectionOperator()->getMapping(),
				*m_solver->getAdvectionOperator()->getDoFHandler(),
				numericVelocity.at(2), ana_uz, local_errors,
				*m_solver->getAdvectionOperator()->getQuadrature(),
				dealii::VectorTools::L2_norm);
		one_component_local_error = local_errors.l2_norm();
		total_local_error += one_component_local_error
				* one_component_local_error;
	}
	const double total_global_error = sqrt(
			dealii::Utilities::MPI::sum(total_local_error,
			MPI_COMM_WORLD));

	m_maxVelocityError = m_solver->m_analyticVelocity.at(0).linfty_norm();
	m_l2VelocityError = total_global_error;

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
