/**
 * @file AdvectionOperator.h
 * @short Abstract class for spatial part of the Advection Operator e_i * dx_i f.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef ADVECTIONOPERATOR_H_
#define ADVECTIONOPERATOR_H_

#include "../utilities/BasicNames.h"
#include "../timeintegration/TimeIntegrator.h"
#include "../solver/DistributionFunctions.h"
#include "../utilities/ConfigNames.h"

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/base/index_set.h"
#include "deal.II/lac/trilinos_sparsity_pattern.h"


namespace natrium {

template<size_t dim>
class ProblemDescription;

class Stencil;

/** @short Abstract class for discretizing the linear advection equation.
 * So far, NATriuM supports a spectral element discontinuous Galerkin solver (SEDGMinLee),
 * and an interpolation-based solver (SemiLagrangian).
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class AdvectionOperator {
protected:

	// Problem Instance
	ProblemDescription<dim>&  m_problem;

	/// Mapping from real space to unit cell
	boost::shared_ptr<dealii::Mapping<dim> > m_mapping;

	/// System matrix L = M^(-1)*(-D+R)
	distributed_sparse_block_matrix m_systemMatrix;

	/// the DQ model (e.g. D2Q9)
	boost::shared_ptr<Stencil> m_stencil;

	/// order of the finite element functions
	size_t m_orderOfFiniteElement;

	/// integration on gauss lobatto nodes
	boost::shared_ptr<dealii::Quadrature<dim> > m_quadrature;

	/// integration on boundary (with gauß lobatto nodes)
	boost::shared_ptr<dealii::Quadrature<dim - 1> > m_faceQuadrature;

	// quadrature on support points with nonsense-weights
	// only for evaluation of values and gradients at the support points
	// not for integration!!
	boost::shared_ptr<dealii::Quadrature<dim> > m_supportPointEvaluation;

	/// Finite Element function on one cell
	boost::shared_ptr<dealii::FiniteElement<dim> > m_fe;

	/// dealii::DoFHandler to distribute the degrees of freedom over the Mesh
	boost::shared_ptr<dealii::DoFHandler<dim> > m_doFHandler;

	dealii::IndexSet m_locallyOwnedDoFs;
	dealii::IndexSet m_locallyRelevantDoFs;

	double m_deltaT;

public:

	/**
	 * @short Constructor
	 * @param problem The ProblemDescription instance that defines the flow
	 * @param fe_order polynomial order of the finite element discretization
	 * @param quad_name name of the quadrature that is used for spatial integration
	 * @param points_name name of the set of support points, as defined in ConfigNames.h
	 * @param stencil the LBM stencil (e.g. D2Q9)
	 * @param delta_t The computational time step. Note that the CFDSolver calls AdvectionOperator::setDeltaT()
	 * 			to set the time step, so as to appy the CFL number to an already refined grid (upon construction
	 * 			of the AdvectionOperator, the grid is not yet refined)
	 * @param dg Flag that indicates whether the finite element is of discontinuous Galerkin type (FE_DGQArbitraryNodes),
	 * 			 otherwise it is a standard Lagrangian finite element (FE_Q)
	 */
	AdvectionOperator(ProblemDescription<dim>& problem,
			size_t fe_order, QuadratureName quad_name,
			SupportPointsName points_name, boost::shared_ptr<Stencil> stencil,
			double delta_t, bool dg);

	/// destructor
	virtual ~AdvectionOperator() {
		m_doFHandler->clear();
	}
	;

	/***
	 * @short Distributes degrees of freedom on the computational grid and extracts locally owned and locally relevant dofs.
	 */
	void distributeDoFs(){
		m_doFHandler->distribute_dofs(*m_fe);
		//get locally owned and locally relevant dofs
		m_locallyOwnedDoFs = m_doFHandler->locally_owned_dofs();
		dealii::DoFTools::extract_locally_relevant_dofs(*m_doFHandler,
				m_locallyRelevantDoFs);
	}

	/**
	 * @short set time step
	 * @param delta_t time step
	 */
	virtual void setDeltaT(double delta_t) {
		m_deltaT = delta_t;
	}

	/////////////////////////////////
	// PURE VIRTUAL FUNCTIONS ///////
	/////////////////////////////////

	/**
	 * @short (re-)assemble system matrix (and system vector, for SEDG). Purely virtual for this class.
	 */
	virtual void reassemble() = 0;

	/**
	 * @short Setup degrees of freedom and fill sparsity pattern of the system matrix. Purely virtual for this class.
	 */
	virtual void setupDoFs() = 0;

	virtual void applyBoundaryConditions(DistributionFunctions& f_old,
			DistributionFunctions& f, double t) = 0;


	/**
	 * @short Stream the distribution functions. Purely virtual for this class.
	 * @param f_old distribution functions at t
	 * @param f distribution functsion at t+delta_t
	 */
	virtual double stream(DistributionFunctions& f_old,
			DistributionFunctions& f, double t) = 0;

	/**
	 * @short Apply boundary conditions. Purely virtual for this class.
	 * @note While the SEDG boundary conditions are directly implemented in the system matrix and vector,
	 * 			the semi-Lagrangian advection operator needs to do something in this functions.
	 */
	//virtual void applyBoundaryConditions(double t) = 0;

	/**
     * @short set the time integrator for the SEDG streaming step.  Purely virtual for this class.
     * @note This function is empty for the semi-Lagrangian streaming, which does not require time integrators.
    */
	virtual void setTimeIntegrator(
			boost::shared_ptr<
					TimeIntegrator<distributed_sparse_block_matrix,
							distributed_block_vector> > timeIntegrator) = 0;

	/**
	 * @short get system vector. To be removed.
	 */
	virtual const distributed_block_vector& getSystemVector() const = 0;



	/////////////////
	// GETTER ///////
	/////////////////
	/**
	 * @short Returns the degree of freedom handler
	 */
	const boost::shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const {
		return m_doFHandler;
	}

	/**
	 * @short return the system matrix
	 */
	const distributed_sparse_block_matrix& getSystemMatrix() const {
		return m_systemMatrix;
	}

	/**
	 * @short return an object to integrate/interpolate on faces, or extract other information
	 * @param flags dealii::UpdateFlags. Depending on which quantity you want to work on,
	 * 			you can specify which quantities are updated at the cell (gradients, shape
	 * 			function values, normal vectors, ...)
	 */
	boost::shared_ptr<dealii::FEFaceValues<dim> > getFEFaceValues(
			const dealii::UpdateFlags & flags) const {
		return boost::make_shared<dealii::FEFaceValues<dim> >(*m_mapping, *m_fe,
				*m_faceQuadrature, flags);
	}

	/**
	 * @short return an object to integrate/interpolate on cells, or extract other information
	 * @param flags dealii::UpdateFlags. Depending on which quantity you want to work on,
	 * 			you can specify which quantities are updated at the cell (gradients, shape
	 * 			function values, normal vectors, ...)
	 */
	boost::shared_ptr<dealii::FEValues<dim> > getFEValues(
			const dealii::UpdateFlags & flags) const {
		return boost::make_shared<dealii::FEValues<dim> >(*m_mapping, *m_fe,
				*m_quadrature, flags);
	}

	/**
	 * @short return the mapping function that maps from an arbitrary quadrilateral cell to a unit cell [0,1]^dim
	 * @note Depending on whether the isCartesian flag is set in your problem description, this mapping can either
	 * 			be a Cartesian or a bilinear ("Q1") mapping.
	 */
	const dealii::Mapping<dim>& getMapping() const {
		return *m_mapping;
	}

	/**
	 * @short return the quadrature object to integrate over faces.
	 */
	const boost::shared_ptr<dealii::Quadrature<dim - 1> >& getFaceQuadrature() const {
		return m_faceQuadrature;
	}

	/***
	 * @short return the finite element object
	 */
	const boost::shared_ptr<dealii::FiniteElement<dim> >& getFe() const {
		return m_fe;
	}

	/**
	 * @short number of degrees of freedom per cell
	 * @note so far, natrium does not support different finite elements on different cells. Thus, the number of dofs per cell is globally well-defined.
	 */
	size_t getNumberOfDoFsPerCell() const {
		return m_fe->dofs_per_cell;
	}


	/**
	 * @short return the quadrature object to integrate over cells.
	 */
	const boost::shared_ptr<dealii::Quadrature<dim> >& getQuadrature() const {
		return m_quadrature;
	}

	/**
	 * @short return the polynomial order of the finite element shape functions
	 */
	size_t getOrderOfFiniteElement() const {
		return m_orderOfFiniteElement;
	}

	/**
	 * @short return the global number of degrees of freedom
	 */
	size_t getNumberOfDoFs() const {
		return getSystemMatrix().block(0, 0).n();
	}

	/**
	 * @short return the set of locally owned degrees of freedom.
	 * Each degree of freedom, i.e. each grid point, is owned by one and only one MPI process.
	 * There are of course other degrees of freedom that are relevant to the calculations (mainly at boundaries of the
	 * locally owned subdomain), those are denoted as locally relevant degrees of freedom.
	 * See the Deal.II glossary for more information.
	 */
	const dealii::IndexSet& getLocallyOwnedDofs() {
		return m_locallyOwnedDoFs;
	}

	/**
	 * @short return the set of locally relevant degrees of freedom.
	 * Each degree of freedom, i.e. each grid point, is owned by one and only one MPI process.
	 * There are of course other degrees of freedom that are relevant to the calculations (mainly at boundaries of the
	 * locally owned subdomain), those are denoted as locally relevant degrees of freedom.
	 * See the Deal.II glossary for more information.
	 */
	const dealii::IndexSet& getLocallyRelevantDofs() {
		return m_locallyRelevantDoFs;
	}

	/**
	 * @short return the computational grid
	 */
	const boost::shared_ptr<Mesh<dim> >& getMesh() const {
		return m_problem.getMesh();
	}

	/**
	 * @short return the LBM stencil (e.g. D2Q9)
	 */
	const boost::shared_ptr<Stencil>& getStencil() const {
		return m_stencil;
	}

	/**
	 * @short return the boundary collection
	 */
	const boost::shared_ptr<BoundaryCollection<dim> >& getBoundaries() const {
		return m_problem.getBoundaries();
	}

	/**
	 * @short estimate the memory consumption of the sparsity pattern
	 * For SEDG, this is 0, because the sparsity pattern is a part of the trilinos matrix.
	 */
	virtual size_t memory_consumption_sparsity_pattern() const {
		return 0;
	}

	/**
	 * @short map degrees of freedom to support points
	 * @param supportPoints. This map is filled with the dof indices and their respective support points.
	 * 			Note that natrium operates on Lagrangian finite elements; thus, support points do always exist.
	 * @note this function just calls dealii::DoFTools::map_dofs_to_support_points()
	 *
	 */
	void mapDoFsToSupportPoints(
			std::map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
		//assert(supportPoints.size() == this->getNumberOfDoFs());
		dealii::DoFTools::map_dofs_to_support_points(*m_mapping, *m_doFHandler,
				supportPoints);
	}

	/**
	 * @short return the support points in real space on a given cell
	 * @param cell the cell
	 * @points points to be filled by this function. Has to have the right size before calling the function, i.e. the dofs per cell.
	 */
	void getCellSupportPoints(typename dealii::DoFHandler<dim>::cell_iterator& cell,
			vector<dealii::Point<dim> >& points) {
		assert(points.size() == getNumberOfDoFsPerCell());
		const vector<dealii::Point<dim> >& unit_support_points =
				m_fe->get_unit_support_points();
		for (size_t i = 0; i < getNumberOfDoFsPerCell(); i++){
			points.at(i) = m_mapping->transform_unit_to_real_cell(cell, unit_support_points.at(i));
		}
	}

	/**
	 * @short A quadrature instance that is meant for evaluation of the finite element functions at the support points
	 * This quadrature should never be used for integration, only for function evaluations,
	 * as the weights are complete nonsense.
	 */
	const dealii::Quadrature<dim>& getSupportPointEvaluation(){
		return *m_supportPointEvaluation;
	}

	double getDeltaT() const {
		return m_deltaT;
	}

	/*const dealii::IndexSet& getLocallyOwnedDoFs() const {
		return m_locallyOwnedDoFs;
	}

	const dealii::IndexSet& getLocallyRelevantDoFs() const {
		return m_locallyRelevantDoFs;
	}*/

	const ProblemDescription<dim> & getProblem() const {
		return m_problem;
	}
};

} /* namespace natrium */
#endif /* ADVECTIONOPERATOR_H_ */
