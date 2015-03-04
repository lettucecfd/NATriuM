/**
 * @file ExponentialTimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */
#include "ExponentialTimeIntegrator.h"
#include <iostream>
#include <cassert>

#include "../utilities/Logging.h"

namespace natrium {

template<> ExponentialTimeIntegrator<distributed_sparse_matrix,
		distributed_vector>::ExponentialTimeIntegrator(double timeStepSize,
		size_t problemSize) :
		TimeIntegrator<distributed_sparse_matrix, distributed_vector>(
				timeStepSize), identityMatrix(makeIdentityMatrix()), H_m(
				makeMatrix(arnoldiSize, arnoldiSize)), H(
				makeMatrix(arnoldiSize + 2, arnoldiSize + 2)) {

}

template<> ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::ExponentialTimeIntegrator(
		double timeStepSize, size_t problemSize, size_t numberOfBlocks) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize), identityMatrix(makeIdentityMatrix()), H_m(
				makeMatrix(arnoldiSize, arnoldiSize)), H(
				makeMatrix(arnoldiSize + 2, arnoldiSize + 2)) {
}

template<class MATRIX, class VECTOR> double ExponentialTimeIntegrator<MATRIX,
		VECTOR>::step(VECTOR& f, const MATRIX& systemMatrix,
		const VECTOR& systemVector, double t, double dt) {

	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());

	if ((0.0 != dt) and dt != this->getTimeStepSize()) {
		this->setTimeStepSize(dt);
		LOG(BASIC) << "Time step size set to " << dt << endl;
	} else if (0.0 == dt) {
		dt = this->getTimeStepSize();
	}

#ifdef WITH_TRILINOS
	if (m_w.size() != f.size()) {
		m_w.reinit(f);
	}

	if (m_v_j.size() != f.size()) {
		m_v_j.reinit(f);
	}

	if (m_v_i.size() != f.size()) {
		m_v_i.reinit(f);
	}
#else
	if (m_w.size() != f.size()) {
		m_w.reinit(f.size(), true);
	}

	if (m_v_j.size() != f.size()) {
		m_v_j.reinit(f.size(), true);
	}

	if (m_v_i.size() != f.size()) {
		m_v_i.reinit(f.size(), true);
	}
#endif

	if (firstColumn.size() != arnoldiSize + 1) {
		firstColumn.reinit(arnoldiSize + 1, true);
	}

	if (m_f.size() != f.size()) {
		m_f.reinit(f.size(), true);
	}

	m_w = f;
	m_v_j = f;
	m_v_i = f;

	for (int i = 0; i < arnoldiSize + 2; i++) {
		for (int j = 0; j < arnoldiSize + 2; j++) {
			H(i, j) = 0;
		}
	}

	numeric_matrix V(f.size(), arnoldiSize + 1); // Orthonormal basis V_(m+1)

	systemMatrix.vmult(m_w, f); 	// w = A*f + u
	m_w.add(systemVector);		// w = A*f + u

	double beta = m_w.l2_norm();

	for (int i = 0; i < f.size(); i++) 		// Arnoldi algorithm (first step)

			{
		V.set(i, 0, m_w(i) / beta);

	}

	for (int j = 0; j < arnoldiSize; j++)  	// Arnoldi algorithm (second step)
			{

		for (int i = 0; i < f.size(); i++) {
			m_v_j(i) = V(i, j);
		}

		systemMatrix.vmult(m_w, m_v_j);

		for (int i = 0; i <= j; i++) {
			for (int k = 0; k < f.size(); k++) {
				m_v_i(k) = V(k, i);
			}

			H(i, j) = m_w * m_v_i;
			m_v_i *= H(i, j);
			m_w -= m_v_i;

		}

		H(j + 1, j) = m_w.l2_norm();
		if (H(j + 1, j) != 0) // && j<arnoldiSize-1)
				{
			for (int k = 0; k < f.size(); k++) {
				V(k, j + 1) = m_w(k) / H(j + 1, j);

			}
		}
		H(arnoldiSize + 1, arnoldiSize) = 1;

	}

	for (int i = 0; i < arnoldiSize; i++) // Transform H to H_m
		for (int j = 0; j < arnoldiSize; j++) {
			{
				H_m(i, j) = H(i, j);
			}
		}

	numeric_matrix phiOne(arnoldiSize);
	numeric_matrix phiTwo(arnoldiSize);
	numeric_matrix phiExtended(arnoldiSize + 1);

	phiFunction(H_m, phiOne, phiTwo);

	for (int i = 0; i < arnoldiSize; i++) {
		for (int j = 0; j < arnoldiSize; j++) {
			phiExtended(i, j) = phiOne(i, j);
		}
		phiExtended(i, arnoldiSize) = 0;
	}

	for (int j = 0; j < arnoldiSize; j++) {
		phiExtended(arnoldiSize, j) = H(arnoldiSize, arnoldiSize - 1)
				* this->getTimeStepSize();
		phiExtended(arnoldiSize, j) *= phiTwo(arnoldiSize - 1, j);
	}

	phiExtended(arnoldiSize, arnoldiSize) = 1;

	for (int i = 0; i < arnoldiSize + 1; i++) {
		firstColumn(i) = phiExtended(i, 0);
	}
	V.vmult(m_f, firstColumn);

	m_f *= (dt * beta);

	for (int i = 0; i < f.size(); i++)

	{
		f(i) += m_f(i);
	}

	/*V.vmult(m_f,firstColumn);
	 m_f*=beta;
	 f = m_f;
	 */

	return t + dt;

}

template double ExponentialTimeIntegrator<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix,
		const distributed_vector& systemVector, double t, double dt);
template double ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix,
		const distributed_block_vector& systemVector, double t, double dt);
} /* namespace natrium */
