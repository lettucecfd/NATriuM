/*
 * matrixAnalysis.cpp
 *
 *  Created on: 03.07.2014
 *      Author: kraemer
 */

#include "matrixAnalysis.h"

namespace natrium {

// constructor/destructor
template<> matrixAnalysis<2>::matrixAnalysis(shared_ptr<CFDSolver<2> > solver):
		m_solver(solver) {

}
template<> matrixAnalysis<3>::matrixAnalysis(shared_ptr<CFDSolver<3> > solver):
		m_solver(solver) {

}
template<> matrixAnalysis<2>::~matrixAnalysis();
template<> matrixAnalysis<3>::~matrixAnalysis();


template<size_t dim>
void matrixAnalysis<dim>::writeSpectrum() {
	// compute eigenvalues of A
	vector<std::complex<double> > eigenvalues;
	copmuteSpectrum(m_solver->getProblemDescription()->getAdvectionOperator()->getSystemMatrix(), eigenvalues);
	// scale by time step (=> eigenvalues of dt*A)
	for (size_t i = 0; i < eigenvalues.size(); i++){
		eigenvalues.at(i) *= m_solver->getTimeStep();
	}
	// open file
	std::stringstream s;
	if ( m_solver->getConfiguration()->isSwitchOutputOff()){
		s << "/tmp/NATriuM_spectrum.txt";
	} else {
		s << m_solver->getConfiguration()->getOutputDirectory().c_str()
							<< "/spectrum.txt";
	}
	std::ofstream spectrumFile(s.str());
	spectrumFile << "# (complex) eigenvalues of system matrix dt*A" << endl;
	spectrumFile << "# Re  Im" << endl;
	// write eigenvalues
	for (size_t i = 0; i < eigenvalues.size(); i++){
		spectrumFile << eigenvalues.at(i).real() << " " << eigenvalues.at(i).imag() << endl;
	}
}
template void matrixAnalysis<2>::writeSpectrum();
template void matrixAnalysis<3>::writeSpectrum();


template<size_t dim>
static void matrixAnalysis<dim>::computeSpectrum(const distributed_sparse_block_matrix& matrix,
		vector<std::complex<double> > & eigenvalues) {
	// make sure the matrix is quadratic
	assert (matrix.n() == matrix.m());
	eigenvalues.clear();
	// use the LAPACK routine for eigenvalues
	dealii::LAPACKFullMatrix<double> fullSystemMatrix = 0;
	fullSystemMatrix.copy_from(matrix);
	fullSystemMatrix.compute_eigenvalues();
	// write eigenvalues to vector
	eigenvalues.resize(matrix.n());
	for (size_t i = 0; i < fullSystemMatrix.n_cols();	i++) {
		eigenvalues.at(i) = fullSystemMatrix.eigenvalue(i);
	}
}
template static void matrixAnalysis<2>::computeSpectrum(const distributed_sparse_block_matrix& matrix,
		vector<std::complex<double> > & eigenvalues);
template static void matrixAnalysis<3>::computeSpectrum(const distributed_sparse_block_matrix& matrix,
		vector<std::complex<double> > & eigenvalues);

} /* namespace natrium */
