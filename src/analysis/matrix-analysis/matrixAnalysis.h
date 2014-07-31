/*
 * matrixAnalysis.h
 *
 *  Created on: 03.07.2014
 *      Author: kraemer
 */

#ifndef MATRIXANALYSIS_H_
#define MATRIXANALYSIS_H_

//#include "complex.h"

#include "solver/CFDSolver.h"

#include "utilities/BasicNames.h"

namespace natrium {

template<size_t dim>
class matrixAnalysis {
private:
	shared_ptr<CFDSolver<dim> > m_solver;

public:

	matrixAnalysis(shared_ptr<CFDSolver<dim> > solver);

	virtual ~matrixAnalysis();

	void writeSpectrum();

	static void computeSpectrum(const distributed_sparse_block_matrix& matrix,
			vector<std::complex<double> > & eigenvalues) ;

	void writePseudospectrum() {

	}

}
;

} /* namespace natrium */

#endif /* MATRIXANALYSIS_H_ */
