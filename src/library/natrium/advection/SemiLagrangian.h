/**
 * @file SemiLagrangian.h
 * @short Semi-Lagrangian advection operator
 * @date 29.04.16
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SEMILAGRANGIAN_H_
#define SEMILAGRANGIAN_H_

#include <map>
#include <array>

#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/lac/block_sparsity_pattern.h"
#include "deal.II/base/quadrature_lib.h"

#include "AdvectionOperator.h"
#include "SemiLagrangianTools.h"
#include "../boundaries/SemiLagrangianBoundaryHandler.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "../utilities/Timing.h"
#include "../utilities/Logging.h"


namespace natrium {

/* forward declaration */
class Stencil;

/**
 * @short Exception class for SemiLagrangian advection operator
 */
class SemiLagrangianException: public NATriuMException {
private:
	std::string message;
public:
	SemiLagrangianException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	SemiLagrangianException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~SemiLagrangianException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/** @short This class solves the linear advection equations by a semi-Lagrangian scheme
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class SemiLagrangian: public AdvectionOperator<dim> {

private:

	/// Mesh
	boost::shared_ptr<Mesh<dim> > m_mesh;

	/// Boundary Description
	boost::shared_ptr<BoundaryCollection<dim> > m_boundaries;

	/// integration on gauss lobatto nodes
	boost::shared_ptr<dealii::Quadrature<dim> > m_quadrature;

	/// integration on boundary (with gauß lobatto nodes)
	boost::shared_ptr<dealii::Quadrature<dim - 1> > m_faceQuadrature;

	/// Finite Element function on one cell
	boost::shared_ptr<dealii::FiniteElement<dim> > m_fe;

	/// dealii::DoFHandler to distribute the degrees of freedom over the Mesh
	boost::shared_ptr<dealii::DoFHandler<dim> > m_doFHandler;

	/// Sparsity Pattern of the sparse matrix
	std::vector<std::vector<dealii::TrilinosWrappers::SparsityPattern> > m_sparsityPattern;

	/// Mapping from real space to unit cell
	const dealii::MappingQ<dim> m_mapping;

	/// System matrix L = M^(-1)*(-D+R)
	distributed_sparse_block_matrix m_systemMatrix;

	distributed_block_vector m_systemVector;

	/// the DQ model (e.g. D2Q9)
	boost::shared_ptr<Stencil> m_stencil;

	/// a map, which connects degrees of freedom with their respective quadrature nodes
	/// m_celldof_to_q_index.at(i)[j] is the support node index q of the j-th dof at a cell
	std::map<size_t, size_t> m_celldof_to_q_index;

	/// a set of maps, which connect degrees of freedom with their respective quadrature nodes
	/// m_facedof_to_q_index.at(i)[j] is the support node index q of the j-th dof at face i
	vector<std::map<size_t, size_t> > m_facedof_to_q_index;

	/// the transposed map of m_facedof_to_q_index
	vector<std::map<size_t, size_t> > m_q_index_to_facedof;

	/// order of the finite element functions
	size_t m_orderOfFiniteElement;

	// time step size
	double m_deltaT;

	// locally owned degrees of freedom (for MPI parallelization)
	dealii::IndexSet m_locallyOwnedDofs;

	// locally relevant degrees of freedom (i.e. ghost layer cells)
	dealii::IndexSet m_locallyRelevantDofs;

	SemiLagrangianBoundaryHandler<dim> m_boundaryHandler;

	/**
	 * @short update the sparsity pattern of the system matrix // the sparse matrix
	 */
	void fillSparseObject(bool sparsity_pattern = false);

	/**
	 * @short update the sparsity pattern of the system matrix
	 */
	void updateSparsityPattern();

	/**
	 * @short map degrees of freedom to quadrature node indices on a cell
	 * @note called by the constructor to initialize m_dof_to_q_index
	 */
	std::map<size_t, size_t> map_celldofs_to_q_index() const;

	/**
	 * @short map degrees of freedom to quadrature node indices on the faces
	 * @note called by the constructor to initialize m_dof_to_q_index
	 */
	vector<std::map<size_t, size_t> > map_facedofs_to_q_index() const;

	/**
	 * @short map quadrature node indices on the faces to degrees of freedom
	 * @note called by the constructor to initialize m_q_index_to_facedof
	 */
	vector<std::map<size_t, size_t> > map_q_index_to_facedofs() const;

public:

	/// constructor
	/**
	 * @short Constructor
	 * @param[in] triangulation The global mesh.
	 * @param[in] orderOfFiniteElement The number of nodes element and dimension
	 * @param[in] stencil the DQ model
	 * @param[in] delta_t time step size; if delta_t = 0, the sparsity pattern is not updated during construction
	 */
	SemiLagrangian(boost::shared_ptr<Mesh<dim> > triangulation,
			boost::shared_ptr<BoundaryCollection<dim> > boundaries,
			size_t orderOfFiniteElement, boost::shared_ptr<Stencil> stencil,
			double delta_t);

	/// destructor
	virtual ~SemiLagrangian() {
		m_doFHandler->clear();
	}
	;


	/**
	 * @short Determines which face is crossed first, when moving from one point inside the cell to a point outside.
	 * @param[in] cell iterator to the active cell that contains the point p_inside
	 * @param[in] p_inside the point inside the cell
	 * @param[in] p_outside the point outside  the cell
	 * @param[out] p_boundary the point where the boundary is hit
	 * @param[out] lambda the parameter lambda that solves   p_boundary_unit = lambda * p_outside_unit + (1-lambda) * p_inside_unit
	 * @return face_id, if a face is crossed; -1, if no face is crossed (i.e. the second point is inside the cell)
	 * @note lambda is calculated for the unit cell. In general, p_boundary = lambda * p_outside + (1-lambda) * p_inside does not hold
	 * @note The current implementation does not do anything special at corner nodes. It prefers x over y over z faces. This may lead to problems later on.
	 */
	int faceCrossedFirst(
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			const dealii::Point<dim>& p_inside,
			const dealii::Point<dim>& p_outside, dealii::Point<dim>& p_boundary,
			double* lambda, size_t* child_id);



	/// function to (re-)assemble linear system
	virtual void reassemble();

	virtual void setupDoFs();

	virtual double stream(DistributionFunctions& f_old,
				DistributionFunctions& f){
		assert (&f_old != &f);
		f_old = f;
		m_systemMatrix.vmult(f.getFStream(), f_old.getFStream());
		return m_deltaT;
	}

	virtual void applyBoundaryConditions(double t) {
		//m_boundaryHandler.apply();

	}


	/// get global system matrix
	virtual const distributed_sparse_block_matrix& getSystemMatrix() const {
		return m_systemMatrix;
	}

	virtual void mapDoFsToSupportPoints(
			std::map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {
		//assert(supportPoints.size() == this->getNumberOfDoFs());
		dealii::DoFTools::map_dofs_to_support_points(m_mapping, *m_doFHandler,
				supportPoints);
	}

	virtual const boost::shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const {
		return m_doFHandler;
	}

	const std::vector<std::vector<dealii::TrilinosWrappers::SparsityPattern> >& getBlockSparsityPattern() const {
		return m_sparsityPattern;
	}

	virtual boost::shared_ptr<dealii::FEFaceValues<dim> > getFEFaceValues(const dealii::UpdateFlags & flags) const {
		return boost::make_shared<dealii::FEFaceValues<dim> >(m_mapping, *m_fe, *m_faceQuadrature,
				flags);
	}

	virtual boost::shared_ptr<dealii::FEValues<dim> > getFEValues(const dealii::UpdateFlags & flags) const {
		return boost::make_shared<dealii::FEValues<dim> >(m_mapping, *m_fe, *m_quadrature,
				flags);
	}

	const dealii::MappingQ<dim>& getMapping() const {
		return m_mapping;
	}

	virtual const std::map<size_t, size_t>& getCelldofToQIndex() const {
		return m_celldof_to_q_index;
	}

	const vector<std::map<size_t, size_t> >& getFacedofToQIndex() const {
		return m_facedof_to_q_index;
	}

	virtual const boost::shared_ptr<dealii::Quadrature<dim - 1> >& getFaceQuadrature() const {
		return m_faceQuadrature;
	}

	virtual const boost::shared_ptr<dealii::FiniteElement<dim> >& getFe() const {
		return m_fe;
	}

	virtual size_t getNumberOfDoFsPerCell() const {
		return m_fe->dofs_per_cell;
	}

	virtual const boost::shared_ptr<dealii::Quadrature<dim> >& getQuadrature() const {
		return m_quadrature;
	}

	virtual size_t getOrderOfFiniteElement() const {
		return m_orderOfFiniteElement;
	}

	virtual size_t getNumberOfDoFs() const {
		return getSystemMatrix().block(0, 0).n();
	}

	virtual const distributed_block_vector& getSystemVector() const {
		return m_systemVector;
	}

	const dealii::IndexSet& getLocallyOwnedDofs() {
		return m_locallyOwnedDofs;
	}
	const dealii::IndexSet& getLocallyRelevantDofs() {
		return m_locallyRelevantDofs;
	}

	virtual const vector<std::map<size_t, size_t> >& getQIndexToFacedof() const {
		return m_q_index_to_facedof;
	}

	virtual const boost::shared_ptr<Mesh<dim> >& getMesh() const {
		return m_mesh;
	}

	virtual const boost::shared_ptr<Stencil>& getStencil() const {
		return m_stencil;
	}

	double getDeltaT() const {
		return m_deltaT;
	}

	virtual void setDeltaT(double deltaT) {
		m_deltaT = deltaT;
		updateSparsityPattern();
		m_boundaryHandler.setTimeStep(deltaT);
	}

	virtual size_t memory_consumption_sparsity_pattern () const {
		size_t mem = 0;
		for (size_t i = 0; i < m_sparsityPattern.size(); i++){
			for (size_t j = 0; j < m_sparsityPattern.at(i).size(); j++){
				mem += m_sparsityPattern.at(i).at(j).memory_consumption();
			}
		}
		return mem;
	}

	/**
	 * @short get i-th neighbor of a cell, incorporating periodic boundaries
	 * @param cell An iterator pointing to cell
	 * @param i face index (=neighbor index)
	 * @return m_doFHandler->end(), if cell has no i-th neighbor (e.g. at solid boundary)
	 *
	 */
	 typename dealii::DoFHandler<dim>::cell_iterator getNeighbor(
	 		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
	 		size_t i);


	 /**
	  * @short fill the neighborhood list
	  * @param cell cell
	  * @param neighborhood the neighborhood object
	  * @note the neighborhood incorporates the current cell, all its neighbors,
	  * 		 and their respective neighbors; each cell has only pointer to it in the neighborhood.
	  */
	 void getNeighborhood(
	 		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
	 		Neighborhood<dim>& neighborhood, size_t n_shells = 1);

	 /**
	  * @short recursively search a point in neighborhood, until is found
	  * @param p the point you search for
	  * @param cell The start cell of the recursive search
	  * @return A cell that contains the point p. If the point was not found, the cell pointer will point to DoFHandler.end()
	  */
	 typename dealii::DoFHandler<dim>::active_cell_iterator recursivelySearchInNeighborhood(
	 		const dealii::Point<dim>& p,
	 		typename dealii::DoFHandler<dim>::active_cell_iterator& cell);


	virtual void setTimeIntegrator(boost::shared_ptr<TimeIntegrator<distributed_sparse_block_matrix,
			distributed_block_vector> > ) {

	}

	const SemiLagrangianBoundaryHandler<dim>& getBoundaryHandler() const {
		return m_boundaryHandler;
	}

	/**
	 * @short Checks if a cell is already in the neighborhood list
	 * @param cell the cell
	 * @param neighborhood the neighborhood list
	 * @return true, if cell is already in the list
	 */
	bool isCellInNeighborhood(
			const typename dealii::DoFHandler<dim>::cell_accessor& cell,
			const Neighborhood<dim>& neighborhood) {
		for (size_t i = 0; i < neighborhood.size(); i++) {
			if (cell.id() == neighborhood.at(i)->id()) {
				return true;
			}
		}
		return false;
	}


}
;

} /* namespace natrium */

#endif /* SEMILAGRANGIAN_H_ */
