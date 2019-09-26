/*
 * VmultLimiter.h
 *
 *  Created on: 29.08.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_DATAPROCESSORS_VMULTLIMITER_H_
#define LIBRARY_NATRIUM_DATAPROCESSORS_VMULTLIMITER_H_

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * Enforces that the result of a sparse matrix-vector multiplication is between the minimum and maximum of the source values.
 * set fi = aij*fj, if (fi > min_{aij>0} fj) and (fi < max{aij>0} fj);
 *        = min_{aij>0} fj, if (fi < min_{aij>0} fj);
 *        = max_{aij>0} fj, if (fi > max_{aij>0} fj)
 */
class VmultLimiter {
public:
	VmultLimiter() {
	}

	virtual ~VmultLimiter() {
	}

	static void apply(const dealii::TrilinosWrappers::BlockSparseMatrix& matrix,
			dealii::TrilinosWrappers::MPI::BlockVector& target,
			const dealii::TrilinosWrappers::MPI::BlockVector& source, distributed_vector& shockSensor);
	static void apply(const dealii::TrilinosWrappers::SparseMatrix& matrix,
			dealii::TrilinosWrappers::MPI::Vector& target,
			const dealii::TrilinosWrappers::MPI::Vector& source, distributed_vector& shockSensor);

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_DATAPROCESSORS_VMULTLIMITER_H_ */
