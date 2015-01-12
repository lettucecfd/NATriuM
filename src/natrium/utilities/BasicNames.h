/**
 * @file BasicNames.h
 * @short Definition of the basic typedefs and names used in the Code
 * @date 30.08.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BASICNAMES_H_
#define BASICNAMES_H_

#include <vector>
#include <iostream>

#include <math.h>

#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"

#include "deal.II/base/point.h"

#include "deal.II/numerics/vector_tools.h"

#include "deal.II/lac/block_sparse_matrix.h"
#include "deal.II/lac/trilinos_vector.h"
#include "deal.II/lac/trilinos_block_vector.h"
#include "deal.II/lac/trilinos_sparse_matrix.h"
#include "deal.II/lac/trilinos_block_sparse_matrix.h"

namespace natrium {

/// The following names will be used throughout natrium
/// by #includeing BasicNames.h they can are used by default
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::size_t;
using std::string;

using boost::shared_ptr;
using boost::make_shared;

/// vector for numeric operations
typedef dealii::Vector<double> numeric_vector;
typedef dealii::BlockVector<double> block_vector;

/// matrix for numeric operations
typedef dealii::FullMatrix<double> numeric_matrix;

/// sparse matrix
typedef dealii::BlockSparseMatrix<double> sparse_block_matrix;
typedef dealii::SparseMatrix<double> sparse_matrix;

#define WITH_TRILINOS
#ifdef WITH_TRILINOS
/// vectors and matrices which can be distributed over different cores
typedef dealii::TrilinosWrappers::MPI::Vector distributed_vector;
typedef dealii::TrilinosWrappers::MPI::BlockVector distributed_block_vector;

typedef dealii::TrilinosWrappers::SparseMatrix distributed_sparse_matrix;
typedef dealii::TrilinosWrappers::BlockSparseMatrix distributed_sparse_block_matrix;
#else

/// vector which can be distributed over different cores
typedef numeric_vector distributed_vector;
typedef block_vector distributed_block_vector;

/// matrix which can be distributed over different cores
typedef sparse_matrix distributed_sparse_matrix;
typedef sparse_block_matrix distributed_sparse_block_matrix;
#endif /* WITH_PETSC */

} /* namespace natrium */


#endif/*BASICNAMES_H_*/

