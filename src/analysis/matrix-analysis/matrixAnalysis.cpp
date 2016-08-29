/**
 * @file matrixAnalysis.cpp
 * @short Computation and output of spectra/pseudospectra of streaming matrices
 * @date 01.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "matrixAnalysis.h"

#include <fstream>
#include <sstream>

namespace natrium {

// constructor/destructor
template<> matrixAnalysis<2>::matrixAnalysis(
		boost::shared_ptr<CFDSolver<2> > solver) :
		m_solver(solver) {

}
template<> matrixAnalysis<3>::matrixAnalysis(
		boost::shared_ptr<CFDSolver<3> > solver) :
		m_solver(solver) {

}
template<> matrixAnalysis<2>::~matrixAnalysis() {

}
template<> matrixAnalysis<3>::~matrixAnalysis() {

}

template<size_t dim>
void matrixAnalysis<dim>::writeSpectrum(bool scale_by_timestep) {
	// compute eigenvalues of A
	vector<std::complex<double> > eigenvalues;
	computeSpectrum(m_solver->getAdvectionOperator()->getSystemMatrix(),
			eigenvalues);
	// scale by time step (=> eigenvalues of dt*A)
	if (scale_by_timestep) {
		for (size_t i = 0; i < eigenvalues.size(); i++) {
			eigenvalues.at(i) *= m_solver->getTimeStepSize();
		}
	}
	// open file
	std::stringstream s;
	if (m_solver->getConfiguration()->isSwitchOutputOff()) {
		s << "/tmp/NATriuM_spectrum.txt";
	} else {
		s << m_solver->getConfiguration()->getOutputDirectory().c_str()
				<< "/spectrum.txt";
	}
	std::ofstream spectrumFile(s.str());
	spectrumFile << "# (complex) eigenvalues of system matrix dt*A" << endl;
	spectrumFile << "# Re  Im" << endl;
	// write eigenvalues
	for (size_t i = 0; i < eigenvalues.size(); i++) {
		spectrumFile << eigenvalues.at(i).real() << " "
				<< eigenvalues.at(i).imag() << endl;
	}
}
template void matrixAnalysis<2>::writeSpectrum(bool scale_by_timestep);
template void matrixAnalysis<3>::writeSpectrum(bool scale_by_timestep);

template<size_t dim>
double matrixAnalysis<dim>::computeSpectrum(
		const distributed_sparse_block_matrix& matrix,
		vector<std::complex<double> > & eigenvalues, double perturbation) {
	// make sure the matrix is quadratic
	assert(matrix.n() == matrix.m());
	eigenvalues.clear();
	// copy the Matrix into LAPACK full matrix
	dealii::LAPACKFullMatrix<double> full_system_matrix(matrix.n(), matrix.n());
	full_system_matrix = 0;
	full_system_matrix.copy_from(matrix);

	double abs_max_eigenvalue = computeSpectrum(full_system_matrix, eigenvalues,
			perturbation);

	//return max absolute of eigenvalues
	return abs_max_eigenvalue;
}
template double matrixAnalysis<2>::computeSpectrum(
		const distributed_sparse_block_matrix& matrix,
		vector<std::complex<double> > & eigenvalues, double perturbation);
template double matrixAnalysis<3>::computeSpectrum(
		const distributed_sparse_block_matrix& matrix,
		vector<std::complex<double> > & eigenvalues, double perturbation);

template<size_t dim>
double matrixAnalysis<dim>::computeSpectrum(
		dealii::LAPACKFullMatrix<double>& matrix,
		vector<std::complex<double> > & eigenvalues, double perturbation) {
	// make sure the matrix is quadratic
	assert(matrix.n() == matrix.m());
	eigenvalues.clear();
	// perturb the matrix (in order to compute pseudospectra)
	// set seed of random number generator
	std::srand(std::time(0));
	if (perturbation != 0.0) {
		for (size_t i = 0; i < matrix.n_rows(); i++) {
			for (size_t j = 0; j < matrix.n_cols(); j++) {
				// random number in [-perturbation,perturbation]
				double rnd = 2 * perturbation
						* (((double) std::rand() / (RAND_MAX)) - 0.5);
				matrix(i, j) += rnd;
			}
		}
	}
	matrix.compute_eigenvalues();
	// write eigenvalues to vector
	eigenvalues.resize(matrix.n());
	double abs_max_eigenvalue = 0.0;
	for (size_t i = 0; i < matrix.n_cols(); i++) {
		eigenvalues.at(i) = matrix.eigenvalue(i);
		if (abs(eigenvalues.at(i)) > abs_max_eigenvalue) {
			abs_max_eigenvalue = abs(eigenvalues.at(i));
		}
	}

	//return max absolute of eigenvalues
	return abs_max_eigenvalue;
}
template double matrixAnalysis<2>::computeSpectrum(
		dealii::LAPACKFullMatrix<double>& matrix,
		vector<std::complex<double> > & eigenvalues, double perturbation);
template double matrixAnalysis<3>::computeSpectrum(
		dealii::LAPACKFullMatrix<double>& matrix,
		vector<std::complex<double> > & eigenvalues, double perturbation);

template<size_t dim>
void matrixAnalysis<dim>::writePseudospectrum(bool scale_by_timestep,
		size_t numberOfCycles, double perturbation) {

	// compute pseudoeigenvalues of A
	vector<std::complex<double> > pseudoeigenvalues;
	vector<std::complex<double> > eigenvalues;
	for (size_t i = 0; i < numberOfCycles; i++) {
		computeSpectrum(m_solver->getAdvectionOperator()->getSystemMatrix(),
				eigenvalues, perturbation);
		// append new eigenvalues to pseudoeigenvalues
		pseudoeigenvalues.insert(pseudoeigenvalues.end(), eigenvalues.begin(),
				eigenvalues.end());

	}
	if (scale_by_timestep) {
		// scale by time step (=> pseudoeigenvalues of dt*A)
		for (size_t i = 0; i < pseudoeigenvalues.size(); i++) {
			pseudoeigenvalues.at(i) *= m_solver->getTimeStepSize();
		}
	}
	// open file
	std::stringstream s;
	if (m_solver->getConfiguration()->isSwitchOutputOff()) {
		s << "/tmp/NATriuM_pseudospectrum.txt";
	} else {
		s << m_solver->getConfiguration()->getOutputDirectory().c_str()
				<< "/pseudospectrum.txt";
	}
	std::ofstream pseudospectrumFile(s.str());
	pseudospectrumFile << "# (complex) pseudoeigenvalues of system matrix dt*A"
			<< endl;
	pseudospectrumFile << "# Re  Im" << endl;
	// write pseudoeigenvalues
	for (size_t i = 0; i < pseudoeigenvalues.size(); i++) {
		pseudospectrumFile << pseudoeigenvalues.at(i).real() << " "
				<< pseudoeigenvalues.at(i).imag() << endl;
	}
}
template void matrixAnalysis<2>::writePseudospectrum(bool scale_by_timestep,
		size_t numberOfCycles, double perturbation);
template void matrixAnalysis<3>::writePseudospectrum(bool scale_by_timestep,
		size_t numberOfCycles, double perturbation);

} /* namespace natrium */
