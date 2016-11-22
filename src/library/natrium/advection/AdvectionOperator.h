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
#include "../timeintegration/TimeIntegrator.h"
#include "../solver/DistributionFunctions.h"

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/base/index_set.h"


namespace natrium {

/** @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class AdvectionOperator {
private:
/*
	// Problem Instance
	boost::shared_ptr<ProblemDescription<dim> > m_problem;

	// Computational grid
	const Mesh<dim>& m_mesh;

	// Boundary description
	const BoundaryCollection<dim> & m_boundaries;

	/// Mapping from real space to unit cell
	const dealii::MappingQ<dim> m_mapping;

	/// Sparsity Pattern of the sparse matrix
	dealii::BlockSparsityPattern m_sparsityPattern;

	/// System matrix L = M^(-1)*(-D+R)
	distributed_sparse_block_matrix m_systemMatrix;

	/// the DQ model (e.g. D2Q9)
	boost::shared_ptr<Stencil> m_stencil;

	/// order of the finite element functions
	size_t m_orderOfFiniteElement;

	/// dealii::DoFHandler to distribute the degrees of freedom over the Mesh
	boost::shared_ptr<dealii::DoFHandler<dim> > m_doFHandler;

	/// a map, which connects degrees of freedom with their respective quadrature nodes
	/// m_celldof_to_q_index.at(i)[j] is the support node index q of the j-th dof at a cell
	std::map<size_t, size_t> m_celldof_to_q_index;

	/// a set of maps, which connect degrees of freedom with their respective quadrature nodes
	/// m_facedof_to_q_index.at(i)[j] is the support node index q of the j-th dof at face i
	vector<std::map<size_t, size_t> > m_facedof_to_q_index;

	/// the transposed map of m_facedof_to_q_index
	vector<std::map<size_t, size_t> > m_q_index_to_facedof;


	/// Mesh
	boost::shared_ptr<Mesh<dim> > m_mesh;

	/// Boundary Description
	boost::shared_ptr<BoundaryCollection<dim> > m_boundaries;

	/// integration on gauss lobatto nodes
	boost::shared_ptr<dealii::QGaussLobatto<dim> > m_quadrature;

	/// integration on boundary (with gauß lobatto nodes)
	boost::shared_ptr<dealii::QGaussLobatto<dim - 1> > m_faceQuadrature;

	/// Finite Element function on one cell
	boost::shared_ptr<dealii::FE_DGQArbitraryNodes<dim> > m_fe;

*/
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
	virtual double stream(DistributionFunctions& f_old,
				DistributionFunctions& f) = 0;

	virtual void applyBoundaryConditions(double t)  = 0;

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

	virtual boost::shared_ptr<dealii::FEFaceValues<dim,dim> > getFEFaceValues(const dealii::UpdateFlags &) const = 0;

	virtual boost::shared_ptr<dealii::FEValues<dim,dim> > getFEValues(const dealii::UpdateFlags &) const = 0;

	virtual const boost::shared_ptr<dealii::FiniteElement<dim> >& getFe() const = 0;

	virtual size_t getNumberOfDoFsPerCell() const = 0;

	virtual const boost::shared_ptr<dealii::Quadrature<dim> >& getQuadrature() const = 0;

	virtual const boost::shared_ptr<dealii::Quadrature<dim - 1> >& getFaceQuadrature() const = 0;

	virtual const std::map<size_t, size_t>& getCelldofToQIndex() const = 0;

	virtual const vector<std::map<size_t, size_t> >& getQIndexToFacedof() const = 0;

	virtual size_t getOrderOfFiniteElement() const = 0;

	virtual const dealii::IndexSet& getLocallyOwnedDofs() = 0;

	virtual const dealii::IndexSet& getLocallyRelevantDofs() = 0;

	virtual const boost::shared_ptr<Mesh<dim> >& getMesh() const = 0;

	virtual const boost::shared_ptr<Stencil>& getStencil() const  = 0;

	virtual void setDeltaT(double ){

	}

	virtual size_t memory_consumption_sparsity_pattern () const {
		return 0;
	}

	virtual void setTimeIntegrator(boost::shared_ptr<TimeIntegrator<distributed_sparse_block_matrix,
			distributed_block_vector> > timeIntegrator) = 0;

};

} /* namespace natrium */
#endif /* ADVECTIONOPERATOR_H_ */
