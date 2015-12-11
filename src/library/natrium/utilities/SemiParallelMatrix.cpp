/*
 * SemiParallelMatrix.cpp
 *
 *  Created on: 13.11.2015
 *      Author: akraem3m
 */

#include "SemiParallelMatrix.h"

namespace natrium {

template<>
void reinitVector<distributed_block_vector>(distributed_block_vector& dst,
		const distributed_block_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src);
	}
}
template<>
void reinitVector<distributed_vector>(distributed_vector& dst,
		const distributed_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src);
	}
}
template<>
void reinitVector<numeric_vector>(numeric_vector& dst,
		const numeric_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src.size(), true);
	}
}
template<>
void reinitVector<block_vector>(block_vector& dst,
		const block_vector& src) {
	if (dst.size() != src.size()) {
		dst.reinit(src.size(), true);
	}
}

template<class VECTOR>
SemiParallelMatrix<VECTOR>::SemiParallelMatrix() :
		m_M() {
}

template<class VECTOR>
SemiParallelMatrix<VECTOR>::SemiParallelMatrix(const VECTOR& prototype,
		size_t n) {
	reinit(prototype, n);
}

template<class VECTOR>
void SemiParallelMatrix<VECTOR>::reinit(const VECTOR& prototype, size_t n) {
	m_M.clear();
	for (size_t i = 0; i < n; i++) {
		VECTOR M_i;
		reinitVector(M_i, prototype);
		M_i = 0;
		m_M.push_back(M_i);
	}
}

template<class VECTOR>
void SemiParallelMatrix<VECTOR>::set(size_t i, size_t j, double value) {
	m_M.at(j)(i) = value;
}

template<class VECTOR>
double SemiParallelMatrix<VECTOR>::operator ()(size_t i, size_t j) const {
	return m_M.at(j)(i);
}

template<class VECTOR>
void SemiParallelMatrix<VECTOR>::vmult(VECTOR& dst,
		const numeric_vector& src) const {
	size_t n = m_M.size();
	assert (src.size() == n);

	dst = 0;
	for (size_t i = 0; i < n; i++){
		dst.add(src(i), m_M.at(i));
	}
}

template<class VECTOR>
void SemiParallelMatrix<VECTOR>::getColumn(size_t j, VECTOR& dst) const{
	assert (j < m_M.size());
	dst = m_M.at(j);
}

template<class VECTOR>
VECTOR&  SemiParallelMatrix<VECTOR>::column(size_t j){
	return m_M.at(j);
}

template<class VECTOR>
void SemiParallelMatrix<VECTOR>::setColumn(size_t j, const VECTOR& src){
	assert (j < m_M.size());
	m_M.at(j) = src;
}

template class SemiParallelMatrix<numeric_vector> ;
template class SemiParallelMatrix<block_vector> ;
template class SemiParallelMatrix<distributed_vector> ;
template class SemiParallelMatrix<distributed_block_vector> ;


} /* namespace natrium */
