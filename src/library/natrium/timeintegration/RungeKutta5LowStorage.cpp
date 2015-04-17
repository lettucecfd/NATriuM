/**
 * @file RungeKutta5LowStorage.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "RungeKutta5LowStorage.h"

#include <cassert>

#include "../utilities/Logging.h"


namespace natrium {

/*template <class MATRIX, class VECTOR> RungeKutta5LowStorage<MATRIX, VECTOR>::RungeKutta5LowStorage(double timeStepSize,
 size_t problemSize): TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
 problemSize), m_Af(problemSize) {}*/

#ifdef WITH_TRILINOS
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
				timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(numberOfBlocks), m_Af(numberOfBlocks) {
	for (size_t i = 0; i < m_Af.size(); i++){
		m_Af.block(i).reinit(problemSize);
		m_df.block(i).reinit(problemSize);
	}
	m_Af.collect_sizes();
	m_df.collect_sizes();
}
#else
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
#endif

template<class MATRIX, class VECTOR> double RungeKutta5LowStorage<MATRIX, VECTOR>::step(
		VECTOR& f, const MATRIX& systemMatrix, const VECTOR& systemVector, double t, double dt){
	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());

	if ((0.0 != dt) and dt != this->getTimeStepSize()){
		this->setTimeStepSize(dt);
		LOG(BASIC) << "Time step size set to " << dt << endl;
	} else if (0.0 == dt){
		dt = this->getTimeStepSize();
	}

#ifdef WITH_TRILINOS
	if (m_Af.size() != f.size()) {
		m_Af.reinit(f);
		/*.n_blocks());
		for (size_t i = 0; i < f.n_blocks(); i++){
			m_Af.block(i).reinit(f.block(i).size());
		}
		m_Af.collect_sizes();*/
	}
	if (m_df.size() != f.size()) {
		m_df.reinit(f);
		/*.n_blocks());
		for (size_t i = 0; i < f.n_blocks(); i++){
			m_df.block(i).reinit(f.block(i).size());
		}
		m_df.collect_sizes();*/
	}
#else
	if (m_Af.size() != f.size()) {
		m_Af.reinit(f.size());
	}
	if (m_df.size() != f.size()) {
		m_df.reinit(f.size());
	}
#endif
	// According to Carpenter and Kennedy (1994) - p. 3
	// df = a*df + h*Af
	// f = f + B* df
	// make first step manually
	systemMatrix.vmult(m_Af, f);
	m_Af.add(systemVector);
	m_Af *= this->getTimeStepSize();
	m_df = m_Af;
	m_df *= m_b.at(0);
	f += m_df;
	m_df /= m_b.at(0);
	for (size_t i = 1; i < 5; i++) {
		systemMatrix.vmult(m_Af, f);
		m_Af.add(systemVector);
		m_Af *= this->getTimeStepSize();
		m_df *= m_a.at(i);
		m_df += m_Af;
		m_df *= m_b.at(i);
		f += m_df;
		m_df /= m_b.at(i);
	}

	return t+dt;

}
template double RungeKutta5LowStorage<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix, const distributed_vector& systemVector, double t, double dt);
template double RungeKutta5LowStorage<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix, const distributed_block_vector& systemVector, double t, double dt);

} /* namespace natrium */