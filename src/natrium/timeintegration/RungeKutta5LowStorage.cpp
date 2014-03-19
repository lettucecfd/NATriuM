/**
 * @file RungeKutta5LowStorage.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "RungeKutta5LowStorage.h"

#include <cassert>

// enable + operator for filling vectors
#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

namespace natrium {

/*template <class MATRIX, class VECTOR> RungeKutta5LowStorage<MATRIX, VECTOR>::RungeKutta5LowStorage(double timeStepSize,
 size_t problemSize): TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
 problemSize), m_Af(problemSize) {}*/

template<> RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector>::RungeKutta5LowStorage(
		double timeStepSize, size_t problemSize) :
		TimeIntegrator<distributed_sparse_matrix, distributed_vector>(
				timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
				problemSize), m_Af(problemSize) {
}

template<> RungeKutta5LowStorage<distributed_sparse_block_matrix,
		distributed_block_vector>::RungeKutta5LowStorage(double timeStepSize,
		size_t problemSize, size_t numberOfBlocks) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(numberOfBlocks,
				problemSize), m_Af(numberOfBlocks, problemSize) {
}


template<class MATRIX, class VECTOR> void RungeKutta5LowStorage<MATRIX, VECTOR>::step(
		VECTOR& f, const MATRIX& systemMatrix) {
	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());
	if (m_Af.size() != f.size()) {
		m_Af.reinit(f.size());
	}
	if (m_df.size() != f.size()) {
		m_df.reinit(f.size());
	}
	// According to Carpenter and Kennedy (1994) - p. 3
	// df = a*df + h*Af
	// f = f + B* df
	// make first step manually
	systemMatrix.vmult(m_Af, f);
	m_Af *= this->getTimeStepSize();
	m_df = m_Af;
	m_df *= m_b.at(0);
	f += m_df;
	m_df /= m_b.at(0);
	for (size_t i = 1; i < 5; i++) {
		systemMatrix.vmult(m_Af, f);
		m_Af *= this->getTimeStepSize();
		m_df *= m_a.at(i);
		m_df += m_Af;
		m_df *= m_b.at(i);
		f += m_df;
		m_df /= m_b.at(i);
	}

}
template void RungeKutta5LowStorage<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix);
template void RungeKutta5LowStorage<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix);

} /* namespace natrium */
