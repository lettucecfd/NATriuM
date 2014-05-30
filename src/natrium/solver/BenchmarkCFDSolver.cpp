/*
 * BenchmarkCFDSolver.cpp
 *
 *  Created on: 27.05.2014
 *      Author: kraemer
 */

#include <solver/BenchmarkCFDSolver.h>

#include "deal.II/dofs/dof_tools.h"

namespace natrium {

template<size_t dim> BenchmarkCFDSolver<dim>::BenchmarkCFDSolver(
		shared_ptr<SolverConfiguration> configuration,
		shared_ptr<Benchmark<dim> > problemDescription) :
		CFDSolver<dim>(configuration, problemDescription), m_benchmark(
				problemDescription) {
	m_supportPoints.resize(this->getNumberOfDoFs());
	m_analyticDensity.reinit(this->getNumberOfDoFs(), true);
	for (size_t i = 0; i < dim; i++) {
		m_analyticVelocity.push_back(
				distributed_vector(this->getNumberOfDoFs()));
	}
	dealii::DoFTools::map_dofs_to_support_points(
			this->getAdvectionOperator()->getMapping(),
			*(this->getAdvectionOperator()->getDoFHandler()), m_supportPoints);

// File for errors
	if ((not configuration->isSwitchOutputOff())
			and (configuration->getOutputTableInterval()
					< configuration->getNumberOfTimeSteps())) {
		std::stringstream s;
		s << configuration->getOutputDirectory().c_str() << "/errors_table.txt";
		if (this->getIterationStart() > 0) {
			m_errorsTableFile = make_shared<std::fstream>(s.str().c_str(),
					std::fstream::out | std::fstream::app);
		} else {
			m_errorsTableFile = make_shared<std::fstream>(s.str().c_str(),
					std::fstream::out);
			(*m_errorsTableFile)
					<< "#  i      t         max |u_numeric|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
					<< endl;
		}
	}

} /*BenchmarkCFDSolver constructor */
template BenchmarkCFDSolver<2>::BenchmarkCFDSolver(
		shared_ptr<SolverConfiguration> configuration,
		shared_ptr<Benchmark<2> > problemDescription);
template BenchmarkCFDSolver<3>::BenchmarkCFDSolver(
		shared_ptr<SolverConfiguration> configuration,
		shared_ptr<Benchmark<3> > problemDescription);

template<> ErrorNorms BenchmarkCFDSolver<2>::getErrors() {
	// this function must not be called more often than once per iteration
	// therefor check a marker value that is set by this function (see below)

	// get analytic and numeric values
	m_benchmark->getAllAnalyticVelocities(getTime(), m_analyticVelocity, m_supportPoints);
	m_benchmark->getAllAnalyticDensities(getTime(), m_analyticDensity, m_supportPoints);
	const vector<distributed_vector>& numericVelocity = getVelocity();
	const distributed_vector& numericDensity = getDensity();

	//calculate errors
	ErrorNorms result;
	//#  i      t         max |u_numeric|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
	result.maxVelocity = getMaxVelocityNorm();
	m_analyticDensity.add(-1.0, numericDensity);
	result.maxDensityError = m_analyticDensity.linfty_norm();
	result.l2DensityError = m_analyticDensity.l2_norm();
	// substract numeric from analytic velocity
	m_analyticVelocity.at(0).add(-1.0, numericVelocity.at(0));
	m_analyticVelocity.at(1).add(-1.0, numericVelocity.at(1));
	// calculate squares
	m_analyticVelocity.at(0).scale(m_analyticVelocity.at(0));
	m_analyticVelocity.at(1).scale(m_analyticVelocity.at(1));
	// calculate ||error (pointwise)||^2
	m_analyticVelocity.at(0).add(m_analyticVelocity.at(1));
	// calculate || error (pointwise) ||
	for (size_t i = 0; i < getNumberOfDoFs(); i++) {
		m_analyticVelocity.at(0)(i) = sqrt(m_analyticVelocity.at(0)(i));
	}
	result.maxVelocityError = m_analyticVelocity.at(0).linfty_norm();
	result.l2VelocityError = m_analyticVelocity.at(0).l2_norm();

	// set marker value that indicates that this function has already been called
	// for the present data
	m_analyticVelocity.at(1)(0) = 31415926;
	return result;
} /*getErrors*/

template<> ErrorNorms BenchmarkCFDSolver<3>::getErrors() {
	// this function must not be called more often than once per iteration
	// therefor check a marker value that is set by this function (see below)

	// get analytic and numeric values
	m_benchmark->getAllAnalyticVelocities(getTime(), m_analyticVelocity, m_supportPoints);
	m_benchmark->getAllAnalyticDensities(getTime(), m_analyticDensity, m_supportPoints);
	const vector<distributed_vector>& numericVelocity = getVelocity();
	const distributed_vector& numericDensity = getDensity();

	//calculate errors
	ErrorNorms result;
	//#  i      t         max |u_numeric|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
	result.maxVelocity = getMaxVelocityNorm();
	m_analyticDensity.add(-1.0, numericDensity);
	result.maxDensityError = m_analyticDensity.linfty_norm();
	result.l2DensityError = m_analyticDensity.l2_norm();
	// substract numeric from analytic velocity
	m_analyticVelocity.at(0).add(-1.0, numericVelocity.at(0));
	m_analyticVelocity.at(1).add(-1.0, numericVelocity.at(1));
	m_analyticVelocity.at(2).add(-1.0, numericVelocity.at(2));
	// calculate squares
	m_analyticVelocity.at(0).scale(m_analyticVelocity.at(0));
	m_analyticVelocity.at(1).scale(m_analyticVelocity.at(1));
	m_analyticVelocity.at(2).scale(m_analyticVelocity.at(2));
	// calculate ||error (pointwise)||^2
	m_analyticVelocity.at(0).add(m_analyticVelocity.at(1));
	m_analyticVelocity.at(0).add(m_analyticVelocity.at(2));
	// calculate || error (pointwise) ||
	for (size_t i = 0; i < getNumberOfDoFs(); i++) {
		m_analyticVelocity.at(0)(i) = sqrt(m_analyticVelocity.at(0)(i));
	}
	result.maxVelocityError = m_analyticVelocity.at(0).linfty_norm();
	result.l2VelocityError = m_analyticVelocity.at(0).l2_norm();

	// set marker value that indicates that this function has already been called
	// for the present data
	m_analyticVelocity.at(1)(0) = 31415926;
	return result;
}

template<size_t dim>
void BenchmarkCFDSolver<dim>::output(size_t iteration) {
	CFDSolver<dim>::output(iteration);
	ErrorNorms norms = getErrors();
	// output error table
	//#  i      t         max |u_numeric|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
	if (iteration % this->getConfiguration()->getOutputTableInterval() == 0) {
		(*m_errorsTableFile) << this->getIteration() << " " << this->getTime() << " "
				<< norms.maxVelocity << " " << norms.maxVelocityError << " "
				<< norms.maxDensityError << " " << norms.l2VelocityError << " "
				<< norms.l2DensityError << endl;
	}

} /*output*/
template void BenchmarkCFDSolver<2>::output(size_t iteration);
template void BenchmarkCFDSolver<3>::output(size_t iteration);

} /* namespace natrium */
