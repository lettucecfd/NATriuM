/*
 * matrixAnalysis.h
 *
 *  Created on: 03.07.2014
 *      Author: kraemer
 */

#ifndef MATRIXANALYSIS_H_
#define MATRIXANALYSIS_H_

#include "complex.h"
#include "fstream"

#include "deal.II/lac/lapack_full_matrix.h"

#include "timeintegration/RungeKutta5LowStorage.h"
#include "utilities/BasicNames.h"
#include "problemdescription/ProblemDescription.h"
#include "solver/CFDSolver.h"

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
