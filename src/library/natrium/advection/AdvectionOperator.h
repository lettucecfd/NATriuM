/**
 * @file AdvectionOperator.h
 * @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef ADVECTIONOPERATOR_H_
#define ADVECTIONOPERATOR_H_

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class AdvectionOperator {


public:

	/// constructor
	AdvectionOperator() {
	}
	;

	/// destructor
	virtual ~AdvectionOperator() {
	}
	;

	/// function to (re-)assemble linear system
	virtual void reassemble() = 0;

	/// make streaming step
	virtual void stream() = 0;
	// TODO is blas installed with dealii? installing blas will speed up the streaming step

	virtual const distributed_sparse_block_matrix& getSystemMatrix() const = 0;

	virtual const distributed_block_vector& getSystemVector() const = 0;

	virtual const shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const = 0;

	virtual void mapDoFsToSupportPoints(
			vector<dealii::Point<dim> >& supportPoints) const = 0;

	virtual const dealii::MappingQ<dim>& getMapping() const = 0;

	/** @short save matrices and status to files
	 *  @param[in] directory directory to save the matrix files to
	 *  @throws AdvectionSolverException
	 */
	virtual void saveCheckpoint(const string& directory) const = 0;

	virtual size_t getNumberOfDoFs() const = 0;

	virtual const shared_ptr<dealii::FE_DGQArbitraryNodes<dim> >& getFe() const = 0;

	//virtual size_t getNumberOfDoFsPerCell() const = 0;
	virtual const shared_ptr<dealii::QGaussLobatto<dim> >& getQuadrature() const = 0;

	virtual const std::map<size_t, size_t>& getCelldofToQIndex() const = 0;

#ifdef WITH_TRILINOS
	virtual const dealii::IndexSet& getLocallyOwnedDofs() = 0;
	virtual const dealii::IndexSet& getLocallyRelevantDofs() = 0;
#endif
};

} /* namespace natrium */
#endif /* ADVECTIONOPERATOR_H_ */
