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
#endif /* WITH_PETSC */

} /* namespace natrium */


#endif/*BASICNAMES_H_*/

