/**
 * @file DistributionFunctions.cpp
 * @short
 * @date 09.11.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DistributionFunctions.h"

namespace natrium {

DistributionFunctions::DistributionFunctions(
		const vector<distributed_vector>& f) :
		m_Q(f.size()), m_f0(f.at(0)) {
#ifdef WITH_TRILINOS
	m_fStream.reinit(m_Q);
#else
	m_fStream.reinit(m_Q, m_f0.size());
#endif
	for (size_t i = 1; i < m_Q; i++) {
		m_fStream.block(i - 1).reinit(f.at(i));
		// reinit does only change the size but not the content
	}
	m_fStream.collect_sizes();
	for (size_t i = 1; i < m_Q; i++) {
		for (size_t j = 0; j < f.at(i).size(); j++) {
			if (m_fStream.block(i - 1).in_local_range(j)) {
				m_fStream.block(i - 1)(j) = f.at(i)(j);
			}
		}
	}
}

distributed_vector& DistributionFunctions::at(size_t i) {
	assert(m_Q > 0);
	assert(i < m_Q);
	if (i == 0) {
		return m_f0;
	} else {
		return m_fStream.block(i - 1);
	}
}

const distributed_vector& DistributionFunctions::at(size_t i) const {
	assert(m_Q > 0);
	assert(i < m_Q);
	if (i == 0) {
		return m_f0;
	} else {
		return m_fStream.block(i - 1);
	}
}

#ifdef WITH_TRILINOS_MPI
void DistributionFunctions::reinit(size_t Q, const dealii::IndexSet &local,
		const dealii::IndexSet &relevant, const MPI_Comm &communicator) {
	m_Q = Q;
	m_f0.reinit(local, relevant, communicator);
	m_fStream.reinit(Q - 1);
	for (size_t i = 0; i < Q - 1; i++) {
		m_fStream.block(i).reinit(m_f0);
	}
	m_fStream.collect_sizes();
}
void DistributionFunctions::reinit(size_t Q, const dealii::IndexSet &local,
		const MPI_Comm &communicator) {
	m_Q = Q;
	m_f0.reinit(local, communicator);
	m_fStream.reinit(Q - 1);
	for (size_t i = 0; i < Q - 1; i++) {
		m_fStream.block(i).reinit(m_f0);
	}
	m_fStream.collect_sizes();
}

#else
void reinit(size_t Q, size_t size) {
	m_Q = Q;
	m_f0.reinit(size);
#ifdef WITH_TRILINOS
	m_fStream.reinit(Q - 1);
	for (size_t i = 0; i < Q - 1; i++) {
		m_fStream.block(i).reinit(m_f0);
	}
#else
	m_fStream.reinit(Q - 1, size);
#endif
	m_fStream.collect_sizes();
}
#endif

void DistributionFunctions::compress(
		dealii::VectorOperation::values operation) {
	m_f0.compress(operation);
	m_fStream.compress(operation);

}

void DistributionFunctions::operator=(const DistributionFunctions& other){
	m_f0 = other.getF0();
	m_fStream = other.getFStream();
}

} /* namespace natrium */
