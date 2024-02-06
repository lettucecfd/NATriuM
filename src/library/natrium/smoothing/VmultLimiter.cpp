/*
 * VmultLimiter.cpp
 *
 *  Created on: 29.08.2016
 *      Author: akraem3m
 */

#include "VmultLimiter.h"
#include "../utilities/CFDSolverUtilities.h"
#include <cfloat>
namespace natrium {

void VmultLimiter::apply(const dealii::TrilinosWrappers::BlockSparseMatrix& matrix,
		dealii::TrilinosWrappers::MPI::BlockVector& target,
		const dealii::TrilinosWrappers::MPI::BlockVector& source){
    (void)matrix;
    (void)target;
    (void)source;
    /*
	// iterate over all blocks
	for (size_t i = 0; i < matrix.n_block_rows(); i++){
		dealii::TrilinosWrappers::MPI::Vector & t = target.block(i);
		for (size_t j = 0; j < matrix.n_block_cols(); j++){
			const dealii::TrilinosWrappers::MPI::Vector & s = source.block(j);
			apply(matrix.block(i,j), t,s);
		}
	}
*/
}


void VmultLimiter::apply(const dealii::TrilinosWrappers::SparseMatrix& matrix,
		dealii::TrilinosWrappers::MPI::Vector& target,
		const dealii::TrilinosWrappers::MPI::Vector& source){
    (void)matrix;
    (void)target;
    (void)source; /*
	// Trilinos matrix format, iteration copied from dealii::TrilinosWrappers::SparseMatrix::print()
	const Epetra_CrsMatrix& M = matrix.trilinos_matrix();
	double *values;
	int *indices;
	int num_entries;
	int glob_i;
	int glob_j;
	double target_i;
	double max;
	double min;

	dealii::TrilinosWrappers::MPI::Vector non_local_source;
	non_local_source.import_nonlocal_data_for_fe(matrix,source);

	for (int i = 0; i < M.NumMyRows(); ++i){
		M.ExtractMyRowView(i, num_entries, values, indices);
		if (num_entries == 0){
			continue;
		}
		//global row index
		glob_i = M.GRID(i);
		// extract value
		target_i = target(glob_i);
		max = std::numeric_limits<double>::min();
		min = std::numeric_limits<double>::max();
		for (dealii::TrilinosWrappers::types::int_type j = 0; j < num_entries; ++j){
			//global column index
			glob_j = M.GCID(indices[j]);

			//if (not source.in_local_range(glob_j))
			//	goto theveryend;

			if (fabs(values[j]) > 1e-12){
				if (non_local_source.el(glob_j) > max){
					max = non_local_source.el(glob_j);
				}
				if (non_local_source.el(glob_j) < min){
					min = non_local_source.el(glob_j);
				}
				if ((target_i <= max) and (target_i >= min)){
					break;
				}
			}
		}
		// limiter
		if (target_i > max){
			target(glob_i) = max;
		} else if (target_i < min){
			target(glob_i) = min;
		}
	}
*/
}


} /* namespace natrium */
