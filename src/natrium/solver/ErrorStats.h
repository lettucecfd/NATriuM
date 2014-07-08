/*
 * ErrorStats.h
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#ifndef ERRORSTATS_H_
#define ERRORSTATS_H_

#include <fstream>

#include "PhysicalProperties.h"
//#include "BenchmarkCFDSolver.h"-> must not be in there because CFDBenchmarkSolver includes this module

#include "../utilities/BasicNames.h"

namespace natrium {
template<size_t dim> class BenchmarkCFDSolver;

/**
 * @short Container for error statistics and table file; only for use with BenchmarkCFDSolver
 * @note As this class is a friend of BenchmarkCFDSolver, it can access its private attributes
 */
template<size_t dim>
class ErrorStats {
private:

	BenchmarkCFDSolver<dim> * m_solver;
	shared_ptr<std::fstream> m_errorsTableFile;
	std::string m_filename;
	const bool m_outputOff;

	// information per iteration
	size_t m_iterationNumber;
	double m_time;
	double m_maxVelocityError;
	double m_maxDensityError;
	double m_l2VelocityError;
	double m_l2DensityError;
	double m_maxUAnalytic;

public:
	ErrorStats(BenchmarkCFDSolver<dim> * cfdsolver,
			const std::string tableFileName = "") :
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
			m_errorsTableFile = make_shared<std::fstream>(tableFileName,
					std::fstream::out | std::fstream::app);
		} else {
			m_errorsTableFile = make_shared<std::fstream>(tableFileName,
					std::fstream::out);
			printHeaderLine();
		}
	}
	void printHeaderLine() {
		(*m_errorsTableFile)
				<< "#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
				<< endl;
	}
	void update() {
		// this function must not be called more often than once per iteration
		// as the data for the analytic solution is constantly overwritten
		// therefor check a marker value that is set by this function (see below)

		m_iterationNumber = m_solver->getIteration();
		m_time = m_solver->getTime();
		// get analytic and numeric values
		// TODO: only assign once (see. addAnalyticSolutionToOutput)
		m_solver->m_benchmark->getAllAnalyticVelocities(m_solver->getTime(),
				m_solver->m_analyticVelocity, m_solver->m_supportPoints);
		m_solver->m_benchmark->getAllAnalyticDensities(m_solver->getTime(),
				m_solver->m_analyticDensity, m_solver->m_supportPoints);
		const vector<distributed_vector>& numericVelocity =
				m_solver->getVelocity();
		const distributed_vector& numericDensity = m_solver->getDensity();

		//#  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2
		m_solver->m_analyticDensity.add(-1.0, numericDensity);
		m_maxDensityError = m_solver->m_analyticDensity.linfty_norm();
		m_l2DensityError = m_solver->m_analyticDensity.l2_norm();

		// calculate maximum analytic velocity norm
		size_t n_dofs = m_solver->getNumberOfDoFs();
		double norm_square = 0.0;
		for (size_t i = 0; i < n_dofs; i++) {
			norm_square = m_solver->m_analyticVelocity.at(0)(i)
					* m_solver->m_analyticVelocity.at(0)(i)
					+ m_solver->m_analyticVelocity.at(1)(i)
							* m_solver->m_analyticVelocity.at(1)(i);
			if (dim == 3){
				norm_square += m_solver->m_analyticVelocity.at(2)(i)
									* m_solver->m_analyticVelocity.at(2)(i);
			}
			if (norm_square > m_maxUAnalytic){
				m_maxUAnalytic = norm_square;
			}
		}
		m_maxUAnalytic = sqrt(m_maxUAnalytic);

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
		m_solver->m_analyticVelocity.at(0).add(
				m_solver->m_analyticVelocity.at(1));
		if (dim == 3) {
			m_solver->m_analyticVelocity.at(0).add(
					m_solver->m_analyticVelocity.at(2));
		}
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
	} /*update*/

	void printNewLine() {
		if (not isUpToDate()) {
			update();
		}
		(*m_errorsTableFile) << m_iterationNumber << " " << m_time << " "
				<< m_maxUAnalytic << " " << m_maxVelocityError << " "
				<< m_maxDensityError << " " << m_l2VelocityError << " "
				<< m_l2DensityError << endl;
	}

	bool isUpToDate() const {
		return (m_iterationNumber == m_solver->getIteration());
	}

	const shared_ptr<std::fstream>& getErrorsTableFile() const {
		return m_errorsTableFile;
	}

	const std::string& getFilename() const {
		return m_filename;
	}

	size_t getIterationNumber() const {
		return m_iterationNumber;
	}

	double getL2DensityError() const {
		return m_l2DensityError;
	}

	double getL2VelocityError() const {
		return m_l2VelocityError;
	}

	double getMaxDensityError() const {
		return m_maxDensityError;
	}

	double getMaxUAnalytic() const {
		return m_maxUAnalytic;
	}

	double getMaxVelocityError() const {
		return m_maxVelocityError;
	}

	double getTime() const {
		return m_time;
	}
};

}/*namespace natrium */

#endif /* ERRORSTATS_H_ */
