/**
 * @file RungeKutta5LowStorage.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "RungeKutta5LowStorage.h"

#include <cassert>

#include "../utilities/Logging.h"
#include "../utilities/SemiParallelMatrix.h"
#include "../utilities/Timing.h"

namespace natrium {

/*template <class MATRIX, class VECTOR> RungeKutta5LowStorage<MATRIX, VECTOR>::RungeKutta5LowStorage(double timeStepSize,
 size_t problemSize): TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
 problemSize), m_Af(problemSize) {}*/

template<> RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector>::RungeKutta5LowStorage(
		double timeStepSize, const distributed_vector& prototype_vector) :
		TimeIntegrator<distributed_sparse_matrix, distributed_vector>(
				timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
				prototype_vector), m_Af(prototype_vector) {
}
template<> RungeKutta5LowStorage<distributed_sparse_block_matrix,
		distributed_block_vector>::RungeKutta5LowStorage(double timeStepSize,
		const distributed_block_vector& prototype_vector) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
				prototype_vector.n_blocks()), m_Af(prototype_vector.n_blocks()) {
	for (size_t i = 0; i < m_Af.size(); i++) {
		m_Af.block(i).reinit(prototype_vector.block(0));
		m_df.block(i).reinit(prototype_vector.block(0));
	}
	m_Af.collect_sizes();
	m_df.collect_sizes();
}
template<> RungeKutta5LowStorage<sparse_matrix, numeric_vector>::RungeKutta5LowStorage(
		double timeStepSize, const numeric_vector& prototype_vector) :
		TimeIntegrator<sparse_matrix, numeric_vector>(timeStepSize), m_a(
				makeA()), m_b(makeB()), m_c(makeC()), m_df(prototype_vector), m_Af(
				prototype_vector) {
}

template<> RungeKutta5LowStorage<sparse_block_matrix, block_vector>::RungeKutta5LowStorage(
		double timeStepSize, const block_vector& prototype_vector) :
		TimeIntegrator<sparse_block_matrix, block_vector>(timeStepSize), m_a(
				makeA()), m_b(makeB()), m_c(makeC()), m_df(
				prototype_vector.n_blocks(), prototype_vector.block(0).size()), m_Af(
				prototype_vector.n_blocks(), prototype_vector.block(0).size()) {
}

template<class MATRIX, class VECTOR> double RungeKutta5LowStorage<MATRIX, VECTOR>::step(
		VECTOR& f, const MATRIX& systemMatrix, VECTOR& systemVector,
		double t, double dt) {
	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());

	if ((0.0 != dt) and dt != this->getTimeStepSize()) {
		this->setTimeStepSize(dt);
		LOG(BASIC) << "Time step size set to " << dt << endl;
	} else if (0.0 == dt) {
		dt = this->getTimeStepSize();
	}

	reinitVector(m_Af, f);
	reinitVector(m_df, f);

	// According to Carpenter and Kennedy (1994) - p. 3
	// df = a*df + h*Af
	// f = f + B* df
	// make first step manually
	{
		TimerOutput::Scope timer_section(Timing::getTimer(), "vmult");
		systemMatrix.vmult(m_Af, f);
	}
	//this->updateSystemVector();
	m_Af += systemVector;
	m_Af *= this->getTimeStepSize();
	m_df = m_Af;
	m_df *= m_b.at(0);
	f += m_df;
	m_df /= m_b.at(0);
	for (size_t i = 1; i < 5; i++) {
		{
			TimerOutput::Scope timer_section(Timing::getTimer(), "vmult");
			systemMatrix.vmult(m_Af, f);
		}
		//this->updateSystemVector();
		m_Af += systemVector;
		m_Af *= this->getTimeStepSize();
		m_df *= m_a.at(i);
		m_df += m_Af;
		m_df *= m_b.at(i);
		f += m_df;
		m_df /= m_b.at(i);
	}

	return t + dt;

}
template double RungeKutta5LowStorage<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix,
		distributed_vector& systemVector, double t, double dt);
template double RungeKutta5LowStorage<sparse_matrix, numeric_vector>::step(
		numeric_vector& f, const sparse_matrix& systemMatrix,
		numeric_vector& systemVector, double t, double dt);
template double RungeKutta5LowStorage<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, double t, double dt);
template double RungeKutta5LowStorage<sparse_block_matrix, block_vector>::step(
		block_vector& f, const sparse_block_matrix& systemMatrix,
		block_vector& systemVector, double t, double dt);

} /* namespace natrium */
