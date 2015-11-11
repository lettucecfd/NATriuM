/**
 * @file BasicNames.h
 * @short Definition of the basic typedefs and names used in the Code;
 * @note As this file is used by most of the others it can contain different compiler flags and other global settings
 * @date 30.08.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BASICNAMES_H_
#define BASICNAMES_H_

// GLOBAL COMPILER FLAGS
// #define WITH_TRILINOS
#define WITH_TRILINOS_MPI
//TODO(AK) merge WITH_TRILINOS_MPI into WITH_TRILINOS

// WITH_TRILINOS_MPI flag includes WITH_TRILINOS flag
#ifdef WITH_TRILINOS_MPI
#ifndef WITH_TRILINOS
#define WITH_TRILINOS
#endif
#endif

#include <vector>
#include <map>
#include <iostream>

#include <math.h>
#include <mpi.h>

#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"


// other deal.II includes
#include "deal.II/base/point.h"
#include "deal.II/numerics/vector_tools.h"
#include "deal.II/lac/block_sparse_matrix.h"
#include "deal.II/base/types.h"
#include "deal.II/base/conditional_ostream.h"


#ifdef WITH_TRILINOS
#include "deal.II/lac/trilinos_vector.h"
#include "deal.II/lac/trilinos_block_vector.h"
#include "deal.II/lac/trilinos_sparse_matrix.h"
#include "deal.II/lac/trilinos_block_sparse_matrix.h"
#include "deal.II/base/index_set.h"
#endif

#ifdef WITH_TRILINOS_MPI
#include "deal.II/distributed/tria.h"
#else
#include "deal.II/grid/tria.h"
#endif

namespace natrium {

/// Typdef Mesh
#ifdef WITH_TRILINOS_MPI
//alias template; works only in C++11
template<size_t dim> using Mesh = dealii::parallel::distributed::Triangulation<dim>;
#else
//alias template; works only in C++11
template<size_t dim> using Mesh = dealii::Triangulation<dim>;
#endif

/// The following names will be used throughout natrium
/// by #including BasicNames.h they can are used by default
using std::vector;
using std::map;
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

// Matrix and Vector classes
#ifdef WITH_TRILINOS
typedef dealii::TrilinosWrappers::SparseMatrix distributed_sparse_matrix;
typedef dealii::TrilinosWrappers::BlockSparseMatrix distributed_sparse_block_matrix;
#else
/// matrix which can be distributed over different cores
typedef sparse_matrix distributed_sparse_matrix;
typedef sparse_block_matrix distributed_sparse_block_matrix;
#endif

#ifdef WITH_TRILINOS_MPI
/// vectors which can be distributed over multiple cores
typedef dealii::TrilinosWrappers::MPI::Vector distributed_vector;
typedef dealii::TrilinosWrappers::MPI::BlockVector distributed_block_vector;

#else
#ifdef WITH_TRILINOS
typedef dealii::TrilinosWrappers::Vector distributed_vector;
typedef dealii::TrilinosWrappers::BlockVector distributed_block_vector;

#else
/// vector which can't be distributed over different cores
typedef numeric_vector distributed_vector;
typedef block_vector distributed_block_vector;
#endif
#endif

// macro to avoid ifdefs when creating a distributed_vector that is actually not distributed
#ifdef WITH_TRILINOS_MPI
// small function to quickly fill in the vectors
#define UNDISTRIBUTED_VECTOR( name, size ) distributed_vector (name) (dealii::complete_index_set(size), MPI_COMM_WORLD)
#else
#define UNDISTRIBUTED_VECTOR( name, size ) distributed_vector (name) ( size )
#endif
// similar for block vector reinit (but no macros needed)
inline void REINIT_UNDISTRIBUTED_BLOCK_VECTOR(distributed_block_vector& bv,
		size_t n_blocks, size_t block_size) {
#ifdef WITH_TRILINOS_MPI
	bv.reinit(n_blocks);
	for (size_t i = 0; i < n_blocks; i++) {
		bv.block(i).reinit(dealii::complete_index_set(block_size), MPI_COMM_WORLD);
	}
	bv.collect_sizes();
#else
#ifdef WITH_TRILINOS
	bv.reinit(n_blocks);
	for (size_t i = 0; i < n_blocks; i++) {
		bv.block(i).reinit(block_size);
	}
	bv.collect_sizes();
#else
	bv.reinit(n_blocks, block_size);
#endif
#endif
} /* REINIT_UNDISTRIBUTED_BLOCK_VECTOR */


/// Parallel cout and cerr
/// Their activity is changed, when MPIGuard is initialized
/// Defined in MPIGuard.cpp
extern dealii::ConditionalOStream perr;
extern dealii::ConditionalOStream pout;


/// Synchronize all MPI processes
inline void MPI_sync(){

	if (!dealii::Utilities::MPI::job_supports_mpi()){
		return;
	}
	//int i = 0;
	// TODO more efficient barrier.
	// sync all MPI processes (barrier)
	//dealii::Utilities::MPI::min_max_avg(i, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

}

/// check if the MPI process has rank zero
/// returns true, if MPI is not used
inline bool is_MPI_rank_0(){
	if (!dealii::Utilities::MPI::job_supports_mpi()){
		return true;
	}
	return (0 == dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
}

}/* namespace natrium */

#endif/*BASICNAMES_H_*/

