/**
 * @file Stats.h
 * @short Contains SolverStats and ErrorStats classes which are responsible for getting data from the solver and storing it to files
 * @date 26.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef STATS_H_
#define STATS_H_

#include <fstream>

#include "PhysicalProperties.h"
//#include "CFDSolver.h"
//#include "BenchmarkCFDSolver.h"

#include "../utilities/BasicNames.h"

namespace natrium {
template<size_t dim> class CFDSolver;
//template<size_t dim> class BenchmarkCFDSolver;

/**
 * @short Container for solver statistics and table file
 * @note As this class is a friend of CFDSolver, it can access its private attributes
 */
template<size_t dim> class SolverStats {
private:

	CFDSolver<dim> * m_solver;
	shared_ptr<std::fstream> m_tableFile;
	std::string m_filename;
	const bool m_outputOff;

	// information per iteration
	size_t m_iterationNumber;
	double m_time;
	double m_maxU;
	double m_kinE;

public:
	SolverStats(CFDSolver<dim> * cfdsolver,
			const std::string tableFileName = "") :
			m_solver(cfdsolver),
			m_filename(tableFileName),
			m_outputOff(tableFileName == ""){
		// be careful: The solver won't be constructed completely when this constructor is called

		// set information
		m_iterationNumber = 1000000000037;
		m_time = 0;
		m_maxU = 0;
		m_kinE = 0;

		if (not m_outputOff) {
			// create file (if necessary)
			if (m_solver->getIterationStart() > 0) {
				m_tableFile = make_shared<std::fstream>(tableFileName,
						std::fstream::out | std::fstream::app);
			} else {
				m_tableFile = make_shared<std::fstream>(tableFileName,
						std::fstream::out);
				printHeaderLine();
			}
		}
	}
	void printHeaderLine() {
		assert(not m_outputOff);
		(*m_tableFile) << "#  i      t      max |u_numeric|    kinE" << endl;
	}
	void update() {
		if (isUpToDate()){
			return;
		}
		m_iterationNumber = m_solver->getIteration();
		m_time = m_solver->getTime();
		m_maxU = m_solver->getMaxVelocityNorm();
		m_kinE = PhysicalProperties<dim>::kineticEnergy(m_solver->getVelocity(),
				m_solver->getDensity());
	}
	void printNewLine() {
		assert(not m_outputOff);
		if (not isUpToDate()) {
			update();
		}
		(*m_tableFile) << m_iterationNumber << " " << m_time << " " << m_maxU
				<< " " << m_kinE << endl;
	}

	size_t getIterationNumber() const {
		return m_iterationNumber;
	}

	double getKinE() const {
		return m_kinE;
	}

	double getMaxU() const {
		return m_maxU;
	}

	shared_ptr<std::fstream> getTableFile() const {
		return m_tableFile;
	}

	double getTime() const {
		return m_time;
	}

	bool isUpToDate() const {
		return (m_iterationNumber == m_solver->getIteration());
	}

	const std::string& getFilename() const {
		return m_filename;
	}
};

/**
 * @short Container for error statistics and table file; only for use with BenchmarkCFDSolver
 * @note As this class is a friend of BenchmarkCFDSolver, it can access its private attributes
 */
/*template<size_t dim>
class ErrorStats {
private:

	const BenchmarkCFDSolver<dim>& m_solver;
	shared_ptr<std::fstream> m_tableFile;

	// information per iteration
	size_t m_iterationNumber;
	double m_time;
	double m_maxVelocityError = 0.0;
	double m_maxDensityError = 0.0;
	double m_l2VelocityError = 0.0;
	double m_l2DensityError = 0.0;
	double m_maxVelocity = 0.0;

public:
	ErrorStats(const BenchmarkCFDSolver<dim>& cfdsolver,
			const std::string tableFileName = "/tmp/natrium_error_table.txt") :
			m_solver(cfdsolver) {
		// set information
		size_t m_iterationNumber = 100000000037;
		double m_time = 0.0;

		// create file (if necessary)
		if (m_solver.getIterationStart() > 0) {
			m_tableFile = make_shared<std::fstream>(tableFileName,
					std::fstream::out | std::fstream::app);
		} else {
			m_tableFile = make_shared<std::fstream>(tableFileName,
					std::fstream::out);
			printHeaderLine();
		}
	}
	void printHeaderLine() {
		(*m_tableFile) << "#  i      t      max |u_numeric|    kinE" << endl;
	}
	void update(){

	}

	void printNewLine() {
		if (not isUpToDate()) {
			update();
		}
		(*m_tableFile) << m_iterationNumber << " " << m_time << " " << endl;
	}

	bool isUpToDate() const {
		return (m_iterationNumber == m_solver.getIteration());
	}
};
*/

} /* namespace natrium */

#endif /* STATS_H_ */
