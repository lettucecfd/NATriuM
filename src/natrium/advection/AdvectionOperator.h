/**
 * @file AdvectionOperator.h
 * @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef ADVECTIONOPERATOR_H_
#define ADVECTIONOPERATOR_H_

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template <size_t dim> class AdvectionOperator {

public:

	/// constructor
	AdvectionOperator(){};

	/// destructor
	virtual ~AdvectionOperator(){};

	/// function to (re-)assemble linear system
	virtual void reassemble() = 0;

	/// make streaming step
	virtual void stream() = 0;
	// TODO is blas installed with dealii? installing blas will speed up the streaming step

	virtual const vector<distributed_sparse_matrix>& getSystemMatrix() const = 0;

	virtual const shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const = 0;

	virtual void mapDoFsToSupportPoints(vector<dealii::Point<dim> >& supportPoints) const = 0;

	size_t getNumberOfDoFs() const {
		return getSystemMatrix().at(0).n();
	}
};

} /* namespace natrium */
#endif /* ADVECTIONOPERATOR_H_ */
