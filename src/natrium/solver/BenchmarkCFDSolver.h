/**
 * @file CFDSolver.h
 * @short CFD Solver, specified for Benchmark applications
 * @date 27.05.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BENCHMARKCFDSOLVER_H_
#define BENCHMARKCFDSOLVER_H_

#include <solver/CFDSolver.h>

#include "deal.II/numerics/data_out.h"

#include "../problemdescription/Benchmark.h"
#include "SolverConfiguration.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @ short an object to store the error norms
 */
struct ErrorNorms {
	double maxVelocityError = 0.0;
	double maxDensityError = 0.0;
	double l2VelocityError = 0.0;
	double l2DensityError = 0.0;
	double maxVelocity = 0.0;
};

/**
 * @short a class that overrides the output function of the CFD solver class with comparisons to a reference solution
 */
template<size_t dim> class BenchmarkCFDSolver: public CFDSolver<dim> {
private:
	/// the problem description, pointed to explicitly as Benchmark
	shared_ptr<Benchmark<dim> > m_benchmark;

	/// support Point of DoFs, in the same order as the DoFs
	vector<dealii::Point<dim> > m_supportPoints;

	/// in order to save allocation time, allocate memory for analytic solution and errors
	distributed_vector m_analyticDensity;
	vector<distributed_vector> m_analyticVelocity;

	/// table out
	shared_ptr<std::fstream> m_errorsTableFile;

public:
	/// constructor
	BenchmarkCFDSolver(shared_ptr<SolverConfiguration> configuration,
			shared_ptr<Benchmark<dim> > problemDescription);

	/**
	 * @short get norms of errors
	 * @note this function calculates the errors in place, which means that it overrides m_analyticDensity/Velocity with garbage
	 */
	ErrorNorms getErrors();

	/**
	 * @short create output data and write to file
	 */
	virtual void output(size_t iteration);

/// gives the possibility for Benchmark instances to add the analytic solution to output
	virtual void addAnalyticSolutionToOutput(dealii::DataOut<dim>& data_out) {
		data_out.add_data_vector(m_analyticDensity, "rho_analytic");
		if (dim == 2) {
			data_out.add_data_vector(m_analyticVelocity.at(0), "ux_analytic");
			data_out.add_data_vector(m_analyticVelocity.at(1), "uy_analytic");
		} else { // dim == 3
			data_out.add_data_vector(m_analyticVelocity.at(0), "ux_analytic");
			data_out.add_data_vector(m_analyticVelocity.at(1), "uy_analytic");
			data_out.add_data_vector(m_analyticVelocity.at(2), "uz_analytic");
		}
	}

/// destructor
	virtual ~BenchmarkCFDSolver() {

	}

	const distributed_vector& getAnalyticDensity() const {
		// check marker value
		if (m_analyticVelocity.at(1)(0) == 31415926){
			throw ConfigurationException("You have to call getAnalyticDensity() before calling getErrors() in a iteration. getErrors() corrupts the analytic density.");
		}
		return m_analyticDensity;
	}

	const vector<distributed_vector>& getAnalyticVelocity() const {
		// check marker value
		if (m_analyticVelocity.at(1)(0) == 31415926){
			throw ConfigurationException("You have to call getAnalyticVelocity() before calling getErrors() in a iteration. getErrors() corrupts the analytic velocity.");
		}
		return m_analyticVelocity;
	}

	const shared_ptr<Benchmark<dim> >& getBenchmark() const {
		return m_benchmark;
	}

	const vector<dealii::Point<2> >& getSupportPoints() const {
		return m_supportPoints;
	}
}
;


} /* namespace natrium */

#endif /* BENCHMARKCFDSOLVER_H_ */

