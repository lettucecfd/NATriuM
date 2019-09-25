/*
 * VmultLimiter.cpp
 *
 *  Created on: 29.08.2016
 *      Author: akraem3m
 */

#include "VmultLimiter.h"
#include <cfloat>
namespace natrium {

void VmultLimiter::apply(const dealii::TrilinosWrappers::BlockSparseMatrix& matrix,
		dealii::TrilinosWrappers::MPI::BlockVector& target,
		const dealii::TrilinosWrappers::MPI::BlockVector& source){
	// iterate over all blocks
	for (size_t i = 0; i < matrix.n_block_rows(); i++){
		dealii::TrilinosWrappers::MPI::Vector & t = target.block(i);
		for (size_t j = 0; j < matrix.n_block_cols(); j++){
			const dealii::TrilinosWrappers::MPI::Vector & s = source.block(j);
			apply(matrix.block(i,j), t,s);
		}
	}

}

void VmultLimiter::apply(const dealii::TrilinosWrappers::SparseMatrix& matrix,
		dealii::TrilinosWrappers::MPI::Vector& target,
		const dealii::TrilinosWrappers::MPI::Vector& source){
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
    dealii::TrilinosWrappers::MPI::Vector local_one(source);

    local_one.import_nonlocal_data_for_fe(matrix,source);
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
            //if (not local_one.in_local_range(glob_j))
             //   cout << "Problem with" << glob_j << " num_entries: " << num_entries << endl;
            //	goto theveryend;
			if (fabs(values[j]) > 1e-12){
                if (local_one(glob_j) > max){
                    max = local_one(glob_j);
				}
                if (local_one(glob_j) < min){
                    min = local_one(glob_j);
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

		theveryend:;

	}

}


} /* namespace natrium */
