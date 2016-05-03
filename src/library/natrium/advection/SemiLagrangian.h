/**
 * @file SemiLagrangian.h
 * @short Semi-Lagrangian advection operator
 * @date 29.04.16
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SEMILAGRANGIAN_H_
#define SEMILAGRANGIAN_H_

#include <map>

#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include <deal.II/dofs/dof_tools.h>
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/lac/block_sparsity_pattern.h"
#include "deal.II/base/quadrature_lib.h"

#include "AdvectionOperator.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"

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
	boost::shared_ptr<dealii::QGaussLobatto<dim> > m_quadrature;

	/// integration on boundary (with gauß lobatto nodes)
	boost::shared_ptr<dealii::QGaussLobatto<dim - 1> > m_faceQuadrature;

	/// Finite Element function on one cell
	boost::shared_ptr<dealii::FE_DGQArbitraryNodes<dim> > m_fe;

	/// dealii::DoFHandler to distribute the degrees of freedom over the Mesh
	boost::shared_ptr<dealii::DoFHandler<dim> > m_doFHandler;

	/// Sparsity Pattern of the sparse matrix
	dealii::BlockSparsityPattern m_sparsityPattern;

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

	/// central flux or Lax-Friedrichs flux (default)
	const bool m_useCentralFlux;

	// locally owned degrees of freedom (for MPI parallelization)
	dealii::IndexSet m_locallyOwnedDofs;

	// locally relevant degrees of freedom (i.e. ghost layer cells)
	dealii::IndexSet m_locallyRelevantDofs;

	/**
	 * @short update the sparsity pattern of the system matrix
	 */
	void updateSparsityPattern();

	/**
	 * @short assemble local mass matrix M
	 * @param[out] massMatrix The mass matrix <phi_i, phi_j>
	 *             For SEDG methods the mass Matrix is diagonal.
	 *
	 */
	// TODO SEDG implemenation with fully diagonal mass matrix
	void assembleLocalMassMatrix(const dealii::FEValues<dim>& feValues,
			size_t dofs_per_cell, vector<double> &massMatrix);

	/**
	 * @short assemble the \f$\alpha\f$-th local derivative matrix
	 * @param[in] coordinate < dim; 0 for Dx, 1 for Dy, 2 for Dz (3D)
	 * @param[out] derivativeMatrix The i-th derivative matrix <D_i phi_j, phi_k>
	 */
	void assembleLocalDerivativeMatrices(const dealii::FEValues<dim>& feValues,
			size_t dofs_per_cell,
			vector<dealii::FullMatrix<double> > &derivativeMatrix) const;

	/**
	 * @short assemble local face matrix;
	 * @param[in] alpha < Q; this matrix is dependent on the flow direction
	 * @param[out] faceMatrix The integral over all faces, incorporating boundary conditions
	 */
	void assembleAndDistributeLocalFaceMatrices(size_t alpha,
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			dealii::FEFaceValues<dim>& feFaceValues,
			dealii::FESubfaceValues<dim>& feSubfaceValues,
			dealii::FEFaceValues<dim>& feNeighborFaceValues,
			const vector<double>& inverseLocalMassMatrix);

	/**
	 * @short calculate system diagonal block matrix  (Dx*eix + Dy*eiy)
	 */
	void calculateAndDistributeLocalStiffnessMatrix(size_t alpha,
			const vector<dealii::FullMatrix<double> > &derivativeMatrices,
			dealii::FullMatrix<double> &systemMatrix,
			const vector<double>& inverseLocalMassMatrix,
			const std::vector<dealii::types::global_dof_index>& globalDoFs,
			size_t dofsPerCell);

	/**
	 * @short assemble and distribute internal face
	 * @param[in] alpha < Q;  the index of the particle flow direction
	 * @param cell the cell to which the face belongs
	 * @param faceNumber the local number of the face (0,1,2 or 3)
	 * @param neighborCell the cell on the other side of the face
	 * @param neighborFaceNumber the local number of the face (0,1,2 or 3), as seen from the neigbor cell
	 * @param feFaceValues is passed in order to avoid allocating memory for each call of the function
	 * @param feSubfaceValues is passed in order to avoid allocating memory for each call of the function
	 * @param feNeighborFaceValues is passed in order to avoid allocating memory for each call of the function
	 */
	void assembleAndDistributeInternalFace(size_t alpha,
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			size_t faceNumber,
			typename dealii::DoFHandler<dim>::cell_iterator& neighborCell,
			size_t neighborFaceNumber, dealii::FEFaceValues<dim>& feFaceValues,
			dealii::FESubfaceValues<dim>& feSubfaceValues,
			dealii::FEFaceValues<dim>& feNeighborFaceValues,
			const vector<double>& inverseLocalMassMatrix);

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

	struct DoFInfo {
		DoFInfo(size_t dof, size_t a, size_t b, const dealii::Point<dim>& x,
				bool is_boundary, dealii::Point<dim> x_b) :
				globalDof(dof), alpha(a), beta(b), oldPoint(x), isBoundary(
						is_boundary), boundaryPoint(x_b) {

		}
		DoFInfo(size_t dof, size_t a, size_t b, const dealii::Point<dim>& x) :
				globalDof(dof), alpha(a), beta(b), oldPoint(x), isBoundary(
						false), boundaryPoint() {

		}
		size_t globalDof; // global degree of freedom
		size_t alpha; // direction id
		size_t beta; // direction one (across boundary)
		dealii::Point<dim> oldPoint; // lagrangian point x^(t-dt)
		bool isBoundary; // flags whether the dof needs boundary handling
		dealii::Point<dim> boundaryPoint; // point where the lagrangian path hits the boundary
	};

	/**
	 * @short A list that stores cell-specific information for assembly
	 */
	typedef std::vector<DoFInfo> XList;

	/**
	 * @short List of neighbors
	 */
	typedef std::vector<typename dealii::DoFHandler<dim>::cell_iterator> Neighborhood;

	/// constructor
	/**
	 * @short Constructor
	 * @param[in] triangulation The global mesh.
	 * @param[in] orderOfFiniteElement The number of nodes element and dimension
	 * @param[in] stencil the DQ model
	 */
	SemiLagrangian(boost::shared_ptr<Mesh<dim> > triangulation,
			boost::shared_ptr<BoundaryCollection<dim> > boundaries,
			size_t orderOfFiniteElement, boost::shared_ptr<Stencil> stencil,
			bool useCentralFlux = false);

	/// destructor
	virtual ~SemiLagrangian() {
		m_doFHandler->clear();
	}
	;

	/**
	 * @short Checks if a cell is already in the neighborhood list
	 * @param cell the cell
	 * @param neighborhood the neighborhood list
	 * @return true, if cell is already in the list
	 */
	bool isCellInNeighborhood(
			const typename dealii::DoFHandler<dim>::cell_accessor& cell,
			const Neighborhood& neighborhood) {
		for (size_t i = 0; i < neighborhood.size(); i++) {
			if (cell.id() == neighborhood.at(i)->id()) {
				return true;
			}
		}
		return false;
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
			size_t i) {
		// cell at periodic boundary
		if (cell->face(i)->at_boundary()) {
			size_t boundaryIndicator = cell->face(i)->boundary_id();
			if (m_boundaries->isPeriodic(boundaryIndicator)) {
				const boost::shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
						m_boundaries->getPeriodicBoundary(boundaryIndicator);
				assert(periodicBoundary->isFaceInBoundary(cell, i));
				typename dealii::DoFHandler<dim>::cell_iterator neighborCell;
				periodicBoundary->getOppositeCellAtPeriodicBoundary(cell,
						neighborCell);
				return neighborCell;
			} else {
				return cell->get_dof_handler().end();
			}
		}
		// cell not at periodic boundary
		return cell->neighbor(i);
	}

	/**
	 * @short fill the neighborhood list
	 * @param cell cell
	 * @param neighborhood the neighborhood object
	 * @note the neighborhood incorporates the current cell, all its neighbors,
	 * 		 and their respective neighbors; each cell has only pointer to it in the neighborhood.
	 */
	void getNeighborhood(
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			Neighborhood& neighborhood) {
		neighborhood.clear();
		// add self
		neighborhood.push_back(cell);
		// as the neighborhood gets extended every time a new neighbor is added
		// this loop runs over all cells (therefor we have to have a break)
		size_t visited = 0;
		size_t n_first_dim_shells = 5;
		if (dim == 3){
			n_first_dim_shells = 25;
		}
		for (size_t c = 0; c < neighborhood.size(); c++) {
			if (visited >= n_first_dim_shells){
				break;
			}
			visited ++;
			for (size_t i = 0; i < dealii::GeometryInfo<dim>::faces_per_cell;
					i++) {
				typename dealii::DoFHandler<dim>::cell_iterator n = getNeighbor(
						neighborhood.at(c), i);
				if (n == cell->get_dof_handler().end()) {
					continue;
				} else if (n->is_artificial()) {
					continue;
				} else if (isCellInNeighborhood(*n, neighborhood)) {
					continue;
				} else if (n->active() or n->is_ghost()) {
					neighborhood.push_back(
							typename dealii::DoFHandler<dim>::active_cell_iterator(
									n));
					continue;
				} else if (n->has_children()) {
					for (size_t j = 0; j < n->n_children(); j++) {
						assert(n->child(j)->active());
						if (isCellInNeighborhood(*(n->child(j)),
								neighborhood)) {
							continue;
						}
						neighborhood.push_back(
								typename dealii::DoFHandler<dim>::active_cell_iterator(
										n->child(j)));
						continue;
					}
				}
			}
		}
	} /* getNeighborhood */

	/**
	 * @short recursively search a point in neighborhood, until is found
	 * @param p the point you search for
	 * @param cell The start cell of the recursive search
	 * @return A cell that contains the point p. If the point was not found, the cell pointer will point to DoFHandler.end()
	 */
	typename dealii::DoFHandler<dim>::active_cell_iterator recursivelySearchInNeighborhood(const dealii::Point<dim>& p,
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell) {
		Neighborhood neighborhood;
		// add self
		neighborhood.push_back(cell);
		for (size_t c = 0; c < neighborhood.size(); c++) {
			if (neighborhood.at(c)->point_inside(p)){
				return neighborhood.at(c);
			}
			for (size_t i = 0; i < dealii::GeometryInfo<dim>::faces_per_cell;
					i++) {
				typename dealii::DoFHandler<dim>::cell_iterator n = getNeighbor(
						neighborhood.at(c), i);
				if (n == cell->get_dof_handler().end()) {
					continue;
				} else if (n->is_artificial()) {
					continue;
				} else if (isCellInNeighborhood(*n, neighborhood)) {
					continue;
				} else if (n->active() or n->is_ghost()) {
					neighborhood.push_back(
							typename dealii::DoFHandler<dim>::active_cell_iterator(
									n));
					continue;
				} else if (n->has_children()) {
					for (size_t j = 0; j < n->n_children(); j++) {
						assert(n->child(j)->active());
						if (isCellInNeighborhood(*(n->child(j)),
								neighborhood)) {
							continue;
						}
						neighborhood.push_back(
								typename dealii::DoFHandler<dim>::active_cell_iterator(
										n->child(j)));
						continue;
					}
				}
			}
		}
		return cell->get_dof_handler().end();
	} /* recursivelySearchInNeighborhood */


	/// function to (re-)assemble linear system
	virtual void reassemble();

	virtual void setupDoFs();

	/// make streaming step
	virtual void stream();

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

	const dealii::BlockSparsityPattern& getBlockSparsityPattern() const {
		return m_sparsityPattern;
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

	const boost::shared_ptr<dealii::QGaussLobatto<dim - 1> >& getFaceQuadrature() const {
		return m_faceQuadrature;
	}

	virtual const boost::shared_ptr<dealii::FE_DGQArbitraryNodes<dim> >& getFe() const {
		return m_fe;
	}

	virtual size_t getNumberOfDoFsPerCell() const {
		return m_fe->dofs_per_cell;
	}

	virtual const boost::shared_ptr<dealii::QGaussLobatto<dim> >& getQuadrature() const {
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
}
;

} /* namespace natrium */

#endif /* SEMILAGRANGIAN_H_ */
