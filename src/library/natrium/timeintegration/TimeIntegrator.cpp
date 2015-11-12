/**
 * @file TimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TimeIntegrator.h"

namespace natrium {

template<class MATRIX, class VECTOR> TimeIntegrator<MATRIX, VECTOR>::TimeIntegrator(
		double timeStepSize) :
		m_timeStepSize(timeStepSize) {
}

template<>
void TimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::reinitVector(distributed_block_vector& dst,
		const distributed_block_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src);
	}
}
template<>
void TimeIntegrator<distributed_sparse_matrix, distributed_vector>::reinitVector(
		distributed_vector& dst, const distributed_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src);
	}
}
template<>
void TimeIntegrator<sparse_matrix, numeric_vector>::reinitVector(
		numeric_vector& dst, const numeric_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src.size(), true);
	}
}
template<>
void TimeIntegrator<sparse_block_matrix, block_vector>::reinitVector(
		block_vector& dst, const block_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src.size(), true);
	}
}


template class TimeIntegrator<distributed_sparse_matrix, distributed_vector>;
template class TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>;
template class TimeIntegrator<sparse_matrix, numeric_vector>;
template class TimeIntegrator<sparse_block_matrix, block_vector>;

/// destructor
//template <class MATRIX, class VECTOR> TimeIntegrator<MATRIX, VECTOR>::~TimeIntegrator(){
//}

} /* namespace natrium */
