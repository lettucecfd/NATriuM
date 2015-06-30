/**
 * @file matrixAnalysis.h
 * @short Header file for the computation of spectra/pseudospectra of streaming matrices
 * @date 01.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef MATRIXANALYSIS_H_
#define MATRIXANALYSIS_H_


#include "natrium/solver/CFDSolver.h"

#include "natrium/utilities/BasicNames.h"

namespace natrium {

/**
 * @short Calculation of spectra and pseudospectra of streaming matrices
 */
template<size_t dim>
class matrixAnalysis {
private:
	shared_ptr<CFDSolver<dim> > m_solver;

public:

	// constructor
	matrixAnalysis(shared_ptr<CFDSolver<dim> > solver);

	// destructor
	virtual ~matrixAnalysis();

	/**
	 * @short Write the spectrum of the CFDSolver's streaming matrix to the file spectrum.txt in the OutputDirectory.
	 * @note If the output directory does not exist: Write to /tmp/NATrium_spectrum.txt
	 * @note The result differs from the spectrum in Min and Lees Paper, because we consider here the full streaming matrix,
	 * which includes all streaming directions.
	 */
	void writeSpectrum();

	/**
	 * @short Write the pseudospectrum of the CFDSolver's streaming matrix to the file pseudospectrum.txt in the OutputDirectory.
	 * The pseudospectrum contains the eigenvalues of all matrices that are "close to the matrix A", i.e. all entries can be
	 * perturbed by a small number.
	 * @param numberOfCycles The number of matrices "close to A", whose eigenvalues will be calculated. Default: 10.
	 * @param perturbation Maximal perturbation of the matrix entries. Default: 0.1.
	 * @note If the output directory does not exist: Write to /tmp/NATrium_pseudospectrum.txt
	 */
	void writePseudospectrum(size_t numberOfCycles = 10, double perturbation = 0.1) ;

	/**
	 * @short Compute the eigenvalues of a sparse block matrix
	 * @param matrix the matrix
	 * @param eigenvalues Vector of eigenvalues. Is automatically resized by the method.
	 * @param perturbation Can be used to perturb all matrix entries by a random number in the interval [-perturbation, perturbation].
	 * This is required for the calculation of pseudospectra. Default: 0.0 (no perturbation).
	 * @note Only practicable for small matrices, as the sparse matrix is turned into a full matrix.
	 * @return maximum of absolute eigenvalues
	 */
	static double computeSpectrum(const distributed_sparse_block_matrix& matrix,
			vector<std::complex<double> > & eigenvalues, double perturbation = 0.0) ;

}
;

} /* namespace natrium */

#endif /* MATRIXANALYSIS_H_ */
