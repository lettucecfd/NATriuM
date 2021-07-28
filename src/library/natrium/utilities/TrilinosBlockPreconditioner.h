/*
 * TrilinosBlockPreconditioner.h
 *
 *  Created on: Feb 12, 2015
 *      Author: kraemer
 */

#ifndef TRILINOSBLOCKPRECONDITIONER_H_
#define TRILINOSBLOCKPRECONDITIONER_H_

#ifdef WITH_TRILINOS

#include "deal.II/base/subscriptor.h"
#include "deal.II/lac/trilinos_parallel_block_vector.h"


namespace natrium {

/// PreconditionIdentity for Block matrices
	class TrilinosBlockPreconditioner: public dealii::Subscriptor {
	public:
		TrilinosBlockPreconditioner() {};
		void vmult (distributed_block_vector &dst,
				const distributed_block_vector &src) const {
			for (size_t i = 0; i < src.n_blocks(); i++){
				dst.block(i) = src.block(i);
				dst.collect_sizes();
			}
		};
		virtual ~TrilinosBlockPreconditioner() {};
	};

} /* namespace natrium */


#endif /* WITH_TRILINOS */

#endif /* TRILINOSBLOCKPRECONDITIONER_H_ */
