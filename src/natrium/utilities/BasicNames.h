/**
 * @file BasicNames.h
 * @short
 * @date 30.08.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef UBLASTYPEDEFS_H_
#define UBLASTYPEDEFS_H_

#include <vector>
#include <iostream>

#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"

#include "deal.II/numerics/vector_tools.h"

#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/lac/petsc_vector.h"
#include "deal.II/lac/petsc_parallel_vector.h"

namespace natrium {

/// The following names will be used throughout natrium
/// by #includeing BasicNames.h they can are used by default
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::size_t;

using boost::shared_ptr;
using boost::make_shared;


/// vector for numeric operations
typedef dealii::Vector<double> numeric_vector;

/// matrix for numeric operations
typedef dealii::FullMatrix<double> numeric_matrix;

/// sparse matrix
typedef dealii::SparseMatrix<double> sparse_matrix;

#undef WITH_PETSC
#ifdef WITH_PETSC
	/// vector which can be distributed over different cores
	typedef dealii::PETScWrappers::MPI::Vector distributed_vector;
	typedef dealii::PETScWrappers::MPI::SparseMatrix distributed_sparse_matrix;
#else
	/// vector which can be distributed over different cores
	typedef numeric_vector distributed_vector;
	typedef sparse_matrix distributed_sparse_matrix;
#endif

/// class which contains basic math functions
class Math {

public:

	/// scalar product
	static double scalar_product(const numeric_vector& x,
			const numeric_vector& y) {
		return x * y;
	}

	/// scale existing vector
	static void scale_vector(double a, numeric_vector& x) {
		x *= a;
	}

	/// scalar times vector
	static numeric_vector scalar_vector(double a, const numeric_vector& x) {
		numeric_vector y(x);
		y *= a;
		return y;
	}

	// add vectors
	static void add_vector(numeric_vector& x, const numeric_vector& y) {
		x += y;
	}

	// subtract vectors
	static void subtract_vector(numeric_vector& x, const numeric_vector& y) {
		x -= y;
	}

	// 2-norm
	static double euclidean_norm(numeric_vector& x) {
		return x.l2_norm();
	}


};

}

/*
 // Copyright (C) 2006 Garth N. Wells.
 // Licensed under the GNU GPL Version 2.
 //
 // Modified by Pawel Krupinski 2007
 // Modified by Andreas Kraemer 2013

 //#define BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
 //#define NDEBUG

 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_of_vector.hpp>
 #include <boost/numeric/ublas/matrix_sparse.hpp>

 // These two files must be included due to a bug in Boost version < 1.33.
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/operation.hpp>

 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

 //#include <boost/numeric/ublas/vector_expression.hpp>

 /// Various typedefs for uBlas data types

 namespace ublas = boost::numeric::ublas;

 // uBlas vector
 typedef ublas::vector<double> ublas_vector;
 typedef ublas::vector_of_vector<double> ublas_vector2d;
 typedef ublas::vector_range<ublas_vector> ublas_vector_range;

 // uBlas dense matrix
 typedef ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> > ublas_dense_matrix_base;
 class ublas_dense_matrix: public ublas_dense_matrix_base {
 public:
 using ublas_dense_matrix_base::operator =;
 ublas_dense_matrix() {
 }
 template<class E>
 ublas_dense_matrix(ublas::matrix_expression<E> const& base) :
 ublas_dense_matrix_base(base) {
 }
 ublas_dense_matrix(const size_t M, const size_t N) :
 ublas_dense_matrix_base(M, N) {
 }
 };
 //typedef ublas_dense_matrix_base ublas_dense_matrix;

 typedef ublas::matrix_range<ublas_dense_matrix_base> ublas_matrix_range;

 // uBlas dense matrix (column major format)
 typedef ublas::matrix<double, ublas::column_major> ublas_matrix_cmajor;
 typedef ublas::matrix_range<ublas_matrix_cmajor> ublas_matrix_cmajor_range;

 // uBlas sparse matrix
 typedef ublas::compressed_matrix<double> ublas_sparse_matrix_base;
 class ublas_sparse_matrix: public ublas_sparse_matrix_base {
 public:
 using ublas_sparse_matrix_base::operator =;
 ublas_sparse_matrix() {
 }
 ublas_sparse_matrix(const size_t M, const size_t N) :
 ublas_sparse_matrix_base(M, N) {
 }
 };

 // uBlas sparse matrix (column major format)
 typedef ublas::compressed_matrix<double, ublas::column_major> ublas_sparse_matrix_cmajor;

 //   // uBlas sparse matrix for temporoary assembly
 //   typedef ublas::generalized_vector_of_vector< double, ublas::row_major,
 //   ublas::vector<ublas::compressed_vector<double> > > ublas_assembly_matrix;
 //
 //   // uBlas sparse matrix for temporoary assembly (column major format)
 //   typedef ublas::generalized_vector_of_vector< double, ublas::column_major,
 //   ublas::vector<ublas::compressed_vector<double> > > ublas_assembly_matrix_cmajor;

 // uBlas upper triangular matrix (column major format)
 typedef ublas::triangular_matrix<double, ublas::upper, ublas::column_major> ublas_matrix_cmajor_tri;
 typedef ublas::matrix_range<ublas_matrix_cmajor_tri> ublas_matrix_cmajor_tri_range;
 typedef ublas::matrix_column<ublas_matrix_cmajor_tri> ublas_matrix_cmajor_tri_column;

 */

#endif /* UBLASTYPEDEFS_H_ */
