/**
 * @file AdvectionOperator.h
 * @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef ADVECTIONOPERATOR_H_
#define ADVECTIONOPERATOR_H_

#include "../utilities/BasicNames.h"
#include "../stencils/Stencil.h"

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/base/index_set.h"


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

	virtual  void setupDoFs() = 0;

	/// make streaming step
	virtual void stream() = 0;
	// TODO is blas installed with dealii? installing blas will speed up the streaming step

	virtual const distributed_sparse_block_matrix& getSystemMatrix() const = 0;

	virtual const distributed_block_vector& getSystemVector() const = 0;

	virtual const boost::shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const = 0;

	virtual void mapDoFsToSupportPoints(
			std::map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const = 0;

	virtual const dealii::MappingQ<dim>& getMapping() const = 0;

	/** @short save matrices and status to files
	 *  @param[in] directory directory to save the matrix files to
	 *  @throws AdvectionSolverException
	 */
	virtual size_t getNumberOfDoFs() const = 0;

	virtual const boost::shared_ptr<dealii::FE_DGQArbitraryNodes<dim> >& getFe() const = 0;

	virtual size_t getNumberOfDoFsPerCell() const = 0;

	virtual const boost::shared_ptr<dealii::QGaussLobatto<dim> >& getQuadrature() const = 0;

	virtual const std::map<size_t, size_t>& getCelldofToQIndex() const = 0;

	virtual const vector<std::map<size_t, size_t> >& getQIndexToFacedof() const = 0;

	virtual size_t getOrderOfFiniteElement() const = 0;

	virtual const dealii::IndexSet& getLocallyOwnedDofs() = 0;

	virtual const dealii::IndexSet& getLocallyRelevantDofs() = 0;

	virtual const boost::shared_ptr<Mesh<dim> >& getMesh() const = 0;

	virtual const boost::shared_ptr<Stencil>& getStencil() const  = 0;

	virtual void setDeltaT(double ){

	}
};

} /* namespace natrium */
#endif /* ADVECTIONOPERATOR_H_ */
