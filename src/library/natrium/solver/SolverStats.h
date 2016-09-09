/*
 * SolverStats.h
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#ifndef SOLVERSTATS_H_
#define SOLVERSTATS_H_

#include <fstream>

#include "PhysicalProperties.h"
//#include "CFDSolver.h" -> must not be in there because CFDSolver includes this module

#include "../utilities/BasicNames.h"

namespace natrium {
template<size_t dim> class CFDSolver;

/**
 * @short Container for solver statistics and table file
 * @note As this class is a friend of CFDSolver, it can access its private attributes
 */
template<size_t dim> class SolverStats {
private:

	CFDSolver<dim> * m_solver;
	boost::shared_ptr<std::fstream> m_tableFile;
	std::string m_filename;
	const bool m_outputOff;

	// information per iteration
	size_t m_iterationNumber;
	double m_time;
	double m_maxU;
	double m_mass;
	double m_kinE;
	double m_enstrophy;
	double m_entropy;
	double m_maxP;
	double m_minP;
	double m_meanVelocityX;

public:
	/**
	 * @short Constructor, if second argument is "" (which it is by default), then no output file is created
	 */
	SolverStats(CFDSolver<dim> * cfdsolver,
			const std::string tableFileName = "") :
			m_solver(cfdsolver), m_filename(tableFileName), m_outputOff(
					tableFileName == "") {
		// be careful: The solver won't be constructed completely when this constructor is called

		// set information
		m_iterationNumber = 1000000000037;
		m_time = 0;
		m_maxU = 0;
		m_mass = 0;
		m_kinE = 0;
		m_enstrophy = 0;
		m_entropy = 0;
		m_maxP = 0;
		m_minP = 0;

		if (not m_outputOff) {
			// create file (if necessary)
			if (m_solver->getIterationStart() > 0) {
				if (is_MPI_rank_0()) {
					m_tableFile = boost::make_shared<std::fstream>(
							tableFileName,
							std::fstream::out | std::fstream::app);
				}
			} else {
				if (is_MPI_rank_0()) {
					m_tableFile = boost::make_shared<std::fstream>(
							tableFileName, std::fstream::out);
				}
				printHeaderLine();
			}
		}

		// sync all MPI processes (barrier)
		MPI_sync();
	}
	void printHeaderLine() {
		assert(not m_outputOff);
		if (is_MPI_rank_0()) {
			(*m_tableFile)
					<< "#  i      t      max |u_numeric|   mass  kinE  enstrophy   entropy   maxP     minP    residuum(rho)   residuum(u)    mean(ux)"
					<< endl;
		}
	}
	void update() {
		if (isUpToDate()) {
			return;
		}
		m_iterationNumber = m_solver->getIteration();
		m_time = m_solver->getTime();
		m_maxU = m_solver->getMaxVelocityNorm();
		m_mass = PhysicalProperties<dim>::mass(m_solver->getDensity(), m_solver->getAdvectionOperator());
		m_kinE = PhysicalProperties<dim>::kineticEnergy(m_solver->getVelocity(),
				m_solver->getDensity(), m_solver->getAdvectionOperator());
		m_enstrophy = PhysicalProperties<dim>::enstrophy(m_solver->getVelocity(), m_solver->getAdvectionOperator());
		m_entropy = PhysicalProperties<dim>::entropy(m_solver->getF(), m_solver->getAdvectionOperator());
		m_maxP = PhysicalProperties<dim>::maximalPressure(
				m_solver->getDensity(),
				m_solver->getStencil()->getSpeedOfSound(), m_minP);
		m_meanVelocityX = PhysicalProperties<dim>::meanVelocityX(
				m_solver->m_velocity.at(0), m_solver->getAdvectionOperator());
		// Residuals must not be calculated because they have to be determined every 10th time steps
	}

	void calulateResiduals(size_t iteration) {
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
					m_solver->m_tmpDensity.linfty_norm());

		}
		// (if first visit of this if-statement, only this part is executed)
		assert(m_solver->m_tmpVelocity.size() == m_solver->m_velocity.size());
		m_solver->m_tmpVelocity.at(0) = m_solver->m_velocity.at(0);
		m_solver->m_tmpVelocity.at(1) = m_solver->m_velocity.at(1);
		m_solver->m_tmpDensity = m_solver->m_density;
	}

	void printNewLine() {
		assert(not m_outputOff);
		if (not isUpToDate()) {
			update();
		}
		if (is_MPI_rank_0()) {
			(*m_tableFile) << m_iterationNumber << " " << m_time << " "
					<< m_maxU << " " << m_mass << " " << m_kinE << " " << m_enstrophy << " " << m_entropy << " " << m_maxP << " " << m_minP
					<< " " << m_solver->m_residuumDensity << " "
					<< m_solver->m_residuumVelocity << " " << m_meanVelocityX
					<< endl;
		}
	}

	size_t getIterationNumber() const {
		return m_iterationNumber;
	}

	double getKinE() const {
		return m_kinE;
	}

	double getMass() const {
		return m_mass;
	}

	double getEnstrophy() const {
		return m_enstrophy;
	}

	double getMaxU() const {
		return m_maxU;
	}

	boost::shared_ptr<std::fstream> getTableFile() const {
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

	double getMaxP() const {
		return m_maxP;
	}

	double getMinP() const {
		return m_minP;
	}

	double getMeanVelocityX() const {
		return m_meanVelocityX;
	}

	double getEntropy() const {
		return m_entropy;
	}

	const bool isOutputOff() const {
		return m_outputOff;
	}
};

} /* namespace natrium */

#endif /* SOLVERSTATS_H_ */
