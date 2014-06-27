/**
 * @file Stats.h
 * @short
 * @date 26.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef STATS_H_
#define STATS_H_

#include "PhysicalProperties.h"
#include "CFDSolver.h"

#include "../utilities/BasicNames.h"

namespace natrium {
template <size_t dim> class CFDSolver;

/**
 * @short Container for solver statistics and table file
 */
template<size_t dim> class SolverStats {
private:

	const CFDSolver<dim>& m_solver;
	shared_ptr<std::fstream> m_tableFile;

	// information per iteration
	size_t m_iterationNumber;
	double m_time;
	double m_maxU;
	double m_kinE;


public:
	SolverStats(const CFDSolver<dim>& cfdsolver, const std::string tableFileName="/tmp/natrium_table.txt"):
		m_solver(cfdsolver){
		// set information
		m_iterationNumber = 1000000000037;
		m_time = 0;
		m_maxU = 0;
		m_kinE = 0;

		// create file (if necessary)
		if (m_solver.getIterationStart() > 0) {
			m_tableFile = make_shared < std::fstream > (tableFileName, std::fstream::out | std::fstream::app);
		} else {
			m_tableFile = make_shared < std::fstream
					> (tableFileName, std::fstream::out);
			printHeaderLine();
		}
	}
	void printHeaderLine() {
		(*m_tableFile) << "#  i      t      max |u_numeric|    kinE" << endl;
	}
	void update(){
		m_iterationNumber = m_solver.getIteration();
		m_time = m_solver.getTime();
		m_maxU = m_solver.getMaxVelocityNorm();
		m_kinE = PhysicalProperties<dim>::kineticEnergy(m_solver.getVelocity(),
						m_solver.getDensity());
	}
	void printNewLine(){
		if (not isUpToDate()){
			update();
		}
		(*m_tableFile) << m_iterationNumber << " " << m_time << " " << m_maxU << " " << m_kinE << endl;
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
		return (m_iterationNumber == m_solver.getIteration());
	}
};

class ErrorStats {

};

} /* namespace natrium */

#endif /* STATS_H_ */
