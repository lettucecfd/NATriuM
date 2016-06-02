/*
 * CompareToBotella.h
 *
 *  Created on: 14.03.2016
 *      Author: akraem3m
 */

#ifndef ANALYSIS_CAVITY_CONVERGENCE_COMPARETOBOTELLA_H_
#define ANALYSIS_CAVITY_CONVERGENCE_COMPARETOBOTELLA_H_

#include "natrium/dataprocessors/DataProcessor.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/utilities/BasicNames.h"


namespace natrium {

class CompareToBotella: public DataProcessor<2> {
private:
	// Reynolds number
	size_t m_Re;
	// y-values
	numeric_vector m_y;
	// u-velocity along vertical line through Geometric center of cavity
	numeric_vector m_referenceU;
	// x-values
	numeric_vector m_x;
	// v-velocity along horizontal line through Geometric center of cavity
	numeric_vector m_referenceV;
	// table filename
	string m_fileName;
	// u_error
	double m_uError;
	double m_vError;
public:
	CompareToBotella(const CFDSolver<2>& solver, size_t reynolds, string filename);
	virtual ~CompareToBotella();
	void makeReferenceU(numeric_vector& u, size_t reynolds_number);
	void makeReferenceV(numeric_vector& v, size_t reynolds_number);
	virtual void apply();
	void printFinalVelocities();

	double getUError() const {
		return m_uError;
	}

	double getVError() const {
		return m_vError;
	}
};

} /* namespace natrium */

#endif /* ANALYSIS_CAVITY_CONVERGENCE_COMPARETOBOTELLA_H_ */
