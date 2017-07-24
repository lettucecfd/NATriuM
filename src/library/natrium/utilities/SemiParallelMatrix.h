/*
 * SemiParallelMatrix.h
 *
 *  Created on: 13.11.2015
 *      Author: akraem3m
 */

#ifndef SEMIPARALLELMATRIX_H_
#define SEMIPARALLELMATRIX_H_

#include "BasicNames.h"


namespace natrium {

/// reinit vector (initialization is different for trilinos data structures,
/// so this functions has two different instantiation)
template<class VECTOR>
void reinitVector(VECTOR& dst, const VECTOR& src);
// forward declare explicit instantiations
template<>
void reinitVector<distributed_block_vector>(distributed_block_vector& dst,
		const distributed_block_vector& src) ;
template<>
void reinitVector<distributed_vector>(distributed_vector& dst,
		const distributed_vector& src);
template<>
void reinitVector<numeric_vector>(numeric_vector& dst,
		const numeric_vector& src) ;
template<>
void reinitVector<block_vector>(block_vector& dst,
		const block_vector& src) ;


/**
 * @short Class that describes a matrix with parallelized columns and non-parallelized rows.
 *
 */
template<class VECTOR> class SemiParallelMatrix {
private:
	vector<VECTOR> m_M;
public:

	/// empty constructor
	SemiParallelMatrix();

	/**
	 * @short Constructor from prototype vector
	 * @param[in] prototype 	A prototype for the matrix columns. Can be dealii::TrilinosWrappers::MPI::Vector
	 * 							distribute the matrix columns over many cores.
	 * @param[in] n				The number of columns
	 */
	SemiParallelMatrix(const VECTOR& prototype, size_t n);

	/**
	 * @short Reinitialize from prototype vector
	 * @param[in] prototype 	A prototype for the matrix columns. Can be dealii::TrilinosWrappers::MPI::Vector
	 * 							distribute the matrix columns over many cores.
	 * @param[in] n				The number of columns
	 */
	void reinit(const VECTOR& prototype, size_t n);

	/**
	 * @short set element (i,j) to value
	 * @param[in] i		Row index
	 * @param[in] j		Column index
	 * @param[in] value Value
	 */
	void set(size_t i, size_t j, double value);

	/**
	 * @short Return matrix element
	 * @param[in] i		Row index
	 * @param[in] j		Column index
	 * @return	  value Value
	 */
	double operator()(size_t i, size_t j) const;

	/**
	 * @short Multiply by vector dst = A*src
	 * @param[out] dst	Result vector. Has the same structure as the columns.
	 * @param[in] src	Vector to be multiplied with. Is not parallelized.
	 */
	void vmult(VECTOR& dst, const numeric_vector& src) const;

	/**
	 * @short extract j-th column of the matrix
	 * @param[in] j 	column index
	 * @param[out] dst	result (column j)
	 *
	 */
	void getColumn(size_t j, VECTOR& dst) const;

	/**
	 * @short get reference to j-th column of the matrix
	 * @param[in] j 		column index
	 * @param[out] return	reference to j-th column
	 *
	 */
	VECTOR&  column(size_t j);

	/**
	 * @short set j-th column of the matrix
	 * @param[in] j 	column index
	 * @param[in] src	new column j
	 *
	 */
	void setColumn(size_t j, const VECTOR& src);
};

} /* namespace natrium */

#endif /*SEMIPARALLELMATRIX_H_ */
