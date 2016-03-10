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
#include "../utilities/Math.h"

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
	boost::shared_ptr<std::fstream> m_errorsTableFile;
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
	double m_l2UAnalytic;

	/**
	 * @short Analytic u function with fixed component, to be passed to integrate_difference
	 */
	template <size_t dim2>
	class AnalyticU: public dealii::Function<dim2>{
	private:
		const dealii::Function<dim2>& m_analyticU;
		size_t m_component;
	public:
		AnalyticU(const dealii::Function<dim2>& ana_u, size_t component):
			m_analyticU(ana_u), m_component(component){
		}
		virtual double value(const dealii::Point<dim2>& x,  const unsigned int component=0) const {
			return m_analyticU.value(x, m_component);
		}
	};

public:
	/**
	 * @short Constructor
	 * @param cfdsolver Instance of CFDSolver object
	 * @param tableFileName Default: "" means: switch output off
	 */
	ErrorStats(BenchmarkCFDSolver<dim> * cfdsolver,
			const std::string tableFileName = "") ;

	/**
	 * @short write header line to table  file
	 */
	void printHeaderLine();

	/**
	 * @short write information of the current iteration to table file
	 */
	void printNewLine();

	/**
	 * @short update errors for the current iteration
	 */
	void update();

	/**
	 * @short check, if errors are up-to-date, i.e. have already been calculated in the current iteration
	 */
	bool isUpToDate() const {
		return (m_iterationNumber == m_solver->getIteration());
	}

	const boost::shared_ptr<std::fstream>& getErrorsTableFile() const {
		return m_errorsTableFile;
	}

	const std::string& getFilename() const {
		return m_filename;
	}

	size_t getIterationNumber() const {
		return m_iterationNumber;
	}

	/**
	 * @short return L2-Error of density,
	 *        \f$ \sqrt{ \sum_{i=1}^{N} (\rho_{i} - \rho_{i}^{ref})^{2} }  \f$
	 * @note The division by the number of dofs is required, because otherwise finer grids result in bigger errors.
	 */
	double getL2DensityError() const {
		return m_l2DensityError;
	}

	/**
	 * @short return L2-Error of velocity,
	 *  \f$ \sqrt{ \sum_{i=1}^{N} \|u_{i} - u_{i}^{ref}\|_{2}^{2} }  \f$
	 * @note The division by the number of dofs is required, because otherwise finer grids result in bigger errors.
	 */
	double getL2VelocityError() const {
		return m_l2VelocityError;
	}

	/**
	 * @short return max error of density
	 *        \f$ max | \rho_{i} - \rho_{i}^{ref} |  \f$
	 */
	double getMaxDensityError() const {
		return m_maxDensityError;
	}


	double getMaxUAnalytic() const {
		return m_maxUAnalytic;
	}

	/**
	 * @short return max error of velocity
	 *        \f$ max  \|u_{i} - u_{i}^{ref}\|_{2}  \f$
	 */
	double getMaxVelocityError() const {
		return m_maxVelocityError;
	}

	double getTime() const {
		return m_time;
	}

	double getL2UAnalytic() const {
		return m_l2UAnalytic;
	}

};

}/*namespace natrium */

#endif /* ERRORSTATS_H_ */
