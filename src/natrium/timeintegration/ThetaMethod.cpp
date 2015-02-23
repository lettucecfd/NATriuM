/**
 * @file ThetaMethod.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "ThetaMethod.h"

#include <cassert>

#include "deal.II/lac/solver_control.h"
#include "deal.II/lac/identity_matrix.h"
#include "deal.II/lac/solver_bicgstab.h"


#ifdef WITH_TRILINOS
#include "deal.II/lac/trilinos_precondition.h"
#include "../utilities/TrilinosBlockPreconditioner.h"
#endif
//#else
#include "deal.II/lac/precondition.h"
//#endif


namespace natrium {

/*template <class MATRIX, class VECTOR> ThetaMethod<MATRIX, VECTOR>::ThetaMethod(double timeStepSize,
 size_t problemSize): TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
 problemSize), m_Af(problemSize) {}*/

template<class MATRIX, class VECTOR> ThetaMethod<MATRIX, VECTOR>::ThetaMethod(
		double timeStepSize, size_t problemSize, double theta) :
		TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_theta(theta) {
}
template ThetaMethod<distributed_sparse_block_matrix, distributed_block_vector>::ThetaMethod(
		double timeStepSize, size_t problemSize, double theta);
template ThetaMethod<distributed_sparse_matrix, distributed_vector>::ThetaMethod(
		double timeStepSize, size_t problemSize, double theta);

#ifdef WITH_TRILINOS
template<> ThetaMethod<distributed_sparse_block_matrix, distributed_block_vector>::ThetaMethod(
		double timeStepSize, size_t problemSize, size_t numberOfBlocks,
		double theta) :
TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
		timeStepSize), m_theta(theta), m_tmpSystemVector(numberOfBlocks) {
	for (size_t i = 0; i < numberOfBlocks; i++) {
		m_tmpSystemVector.block(i).reinit(problemSize);
	}
	m_tmpSystemVector.collect_sizes();
}
#else
template<> ThetaMethod<distributed_sparse_block_matrix, distributed_block_vector>::ThetaMethod(
		double timeStepSize, size_t problemSize, size_t numberOfBlocks,
		double theta) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize), m_theta(theta), m_tmpSystemVector(numberOfBlocks,
				problemSize) {
}
#endif

template<> double ThetaMethod<distributed_sparse_matrix, distributed_vector>::step(
		distributed_vector& f, const distributed_sparse_matrix& systemMatrix,
		const distributed_vector& systemVector, double t, double dt) {
	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());
	assert(systemVector.size() == systemMatrix.n());

#ifdef WITH_TRILINOS
	if (m_tmpSystemVector.size() != f.size()) {
		m_tmpSystemVector.reinit(f);
	}
	// check equality of sparsity patterns
	// check equality of sparsity patterns
	if (m_tmpMatrix.memory_consumption() != systemMatrix.memory_consumption()) {
		m_tmpMatrix.copy_from(systemMatrix);
	}
#else
	if (m_tmpSystemVector.size() != f.size()) {
		m_tmpSystemVector.reinit(f.size());
	}
	if ((m_tmpMatrix.empty()) or
	// the next check should give true if the sparsity patterns are equal
	// and false, else. n_nonzero_elements returns the number of entries
	// in the sparsity pattern, not the actual number of nonzero entries
			(m_tmpMatrix.n_nonzero_elements()
					!= systemMatrix.n_nonzero_elements())) {
		m_tmpMatrix.reinit(systemMatrix.get_sparsity_pattern());
	}
#endif
	// dt*b
	m_tmpSystemVector = systemVector;
	m_tmpSystemVector *= this->getTimeStepSize();
	// dt*A*f(t) + dt*b
	m_tmpMatrix.copy_from(systemMatrix);
	m_tmpMatrix *= this->getTimeStepSize();
	m_tmpMatrix.vmult_add(m_tmpSystemVector, f);
	// I-theta*dt*A
	m_tmpMatrix *= (-m_theta);
	for (size_t i = 0; i < m_tmpMatrix.n(); i++) {
		m_tmpMatrix.add(i, i, 1.0);
	}
	// (I-theta*dt*A)f(t) + dt*A*f(t) + dt*b
	m_tmpMatrix.vmult_add(m_tmpSystemVector, f);

	dealii::SolverControl solver_control(1000, 1e-8, false, false);	//* m_tmpSystemVector.l2_norm());
	dealii::SolverBicgstab<distributed_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
	bicgstab.solve(m_tmpMatrix, f, m_tmpSystemVector,
			dealii::TrilinosWrappers::PreconditionIdentity());
#else
	bicgstab.solve(m_tmpMatrix, f, m_tmpSystemVector,
			dealii::PreconditionIdentity());	//,	           preconditioner);
#endif

}
template<> double ThetaMethod<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix,
		const distributed_block_vector& systemVector, double t, double dt) {
	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());
	assert(systemVector.size() == systemMatrix.n());

#ifdef WITH_TRILINOS
	if (m_tmpSystemVector.size() != f.size()) {
		m_tmpSystemVector.reinit(f);
	}
	// check equality of sparsity patterns
	if (m_tmpMatrix.memory_consumption() != systemMatrix.memory_consumption()) {
		size_t n_blocks = systemMatrix.n_block_rows();
		assert (systemMatrix.n_block_cols() == systemMatrix.n_block_rows());
		m_tmpMatrix.reinit(n_blocks, n_blocks);
		for (size_t I = 0; I < n_blocks; I++) {
			for (size_t J = 0; J < n_blocks; J++) {
				m_tmpMatrix.block(I, J).reinit(systemMatrix.block(I, J));
			}
		}
	}
#else
	if (m_tmpSystemVector.size() != f.size()) {
		m_tmpSystemVector.reinit(f.size());
	}
	if ((m_tmpMatrix.empty()) or
	// the next check should give true if the sparsity patterns are equal
	// and false, else. n_nonzero_elements returns the number of entries
	// in the sparsity pattern, not the actual number of nonzero entries
			(m_tmpMatrix.n_nonzero_elements()
					!= systemMatrix.n_nonzero_elements())) {
		m_tmpMatrix.reinit(systemMatrix.get_sparsity_pattern());
	}
#endif
	// dt*b
	m_tmpSystemVector = systemVector;
	m_tmpSystemVector *= this->getTimeStepSize();
	// dt*A*f(t) + dt*b
	m_tmpMatrix.copy_from(systemMatrix);
	size_t n_blocks = systemMatrix.n_block_rows();
	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			assert(
					m_tmpMatrix.block(I, J).l1_norm()
							== systemMatrix.block(I, J).l1_norm());
		}
	}
	m_tmpMatrix *= this->getTimeStepSize();
	m_tmpMatrix.vmult_add(m_tmpSystemVector, f);
	// I-theta*dt*A
	m_tmpMatrix *= (-m_theta);
	for (size_t I = 0; I < m_tmpMatrix.n_block_cols(); I++){
		for (size_t i = 0; i < m_tmpMatrix.block(I,I).n(); i++) {
			m_tmpMatrix.block(I,I).add(i, i, 1.0);
		}
	}

	// (I-theta*dt*A)f(t) + dt*A*f(t) + dt*b
	m_tmpMatrix.vmult_add(m_tmpSystemVector, f);


	//dealii::PreconditionBlockSSOR<MATRIX> preconditioner(m_tmpMatrix);
	dealii::SolverControl solver_control(1000, 1e-8, false, false);	//* m_tmpSystemVector.l2_norm());
	dealii::SolverBicgstab<distributed_block_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
	bicgstab.solve(m_tmpMatrix, f, m_tmpSystemVector,
			dealii::PreconditionIdentity());//TrilinosBlockPreconditioner());
#else
	bicgstab.solve(m_tmpMatrix, f, m_tmpSystemVector,
			dealii::PreconditionIdentity());	//,	           preconditioner);
#endif

}

} /* namespace natrium */
