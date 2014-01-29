/**
 * @file SEDGMinLee.h
 * @short Advection operator proposed by Min and Lee (2011): A spectral-elemennt discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SEDGMINLEE_H_
#define SEDGMINLEE_H_

#include <map>

#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include <deal.II/dofs/dof_tools.h>
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/base/quadrature_lib.h"

#include "AdvectionOperator.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../boltzmannmodels/BoltzmannModel.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Global data which is used, e.g., by Min and Lee (2011): A spectral-element discontinuous
 *         Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 *         including particle distributions f, system matrix L, diagonal mass matrix M,
 *         gradient matrices Dx, Dy, (Dz) and boundary matrix R
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class SEDGMinLee: public AdvectionOperator<dim> {

private:

	/// Triangulation
	shared_ptr<dealii::Triangulation<dim> > m_tria;

	/// Boundary Description
	shared_ptr<BoundaryCollection<dim> > m_boundaries;

	/// integration on gauss lobatto nodes
	shared_ptr<dealii::QGaussLobatto<dim> > m_quadrature;

	/// integration on boundary (with gauß lobatto nodes)
	shared_ptr<dealii::QGaussLobatto<dim - 1> > m_faceQuadrature;

	/// Finite Element function on one cell
	shared_ptr<dealii::FE_DGQArbitraryNodes<dim> > m_fe;

	/// dealii::DoFHandler to distribute the degrees of freedom over the Triangulation
	shared_ptr<dealii::DoFHandler<dim> > m_doFHandler;

	/// Sparsity Pattern of the sparse matrix
	dealii::SparsityPattern m_sparsityPattern;

	/// Mapping from real space to unit cell
	const dealii::MappingQ1<dim> m_mapping;

	/// global mass matrix M
	distributed_vector m_massMatrix;

	/// System matrix L = M^(-1)*(-D+R)
	vector<distributed_sparse_matrix> m_systemMatrix;

	/// the DQ model (e.g. D2Q9)
	shared_ptr<BoltzmannModel> m_boltzmannModel;

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
			size_t dofs_per_cell, size_t n_q_points, vector<double> &massMatrix,
			const std::vector<dealii::types::global_dof_index>& globalDoFs);

	/**
	 * @short assemble the i-th local derivative matrix
	 * @param[in] coordinate < dim; 0 for Dx, 1 for Dy, 2 for Dz (3D)
	 * @param[out] derivativeMatrix The i-th derivative matrix <D_i phi_j, phi_k>
	 */
	void assembleLocalDerivativeMatrices(const dealii::FEValues<dim>& feValues,
			size_t dofs_per_cell, size_t n_q_points,
			vector<dealii::FullMatrix<double> > &derivativeMatrix) const;

	/**
	 * @short assemble local face matrix;
	 * @param[in] i < Q; this matrix is dependent on the flow direction
	 * @param[out] faceMatrix The integral over all faces, incorporating boundary conditions
	 */
	void assembleAndDistributeLocalFaceMatrices(size_t i,
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			dealii::FEFaceValues<dim>& feFaceValues,
			dealii::FESubfaceValues<dim>& feSubfaceValues,
			dealii::FEFaceValues<dim>& feNeighborFaceValues,
			size_t dofs_per_cell, size_t n_q_points,
			dealii::FullMatrix<double> &faceMatrix);

	/**
	 * @short calculate system diagonal block matrix  (Dx*eix + Dy*eiy)
	 */
	void calculateAndDistributeLocalStiffnessMatrix(size_t i,
			const vector<dealii::FullMatrix<double> > &derivativeMatrices,
			dealii::FullMatrix<double> &systemMatrix,
			const std::vector<dealii::types::global_dof_index>& globalDoFs,
			size_t dofsPerCell);

	/**
	 * @short calculate A <- M^-1 * A
	 * @param[in/out] matrix sparse matrix A
	 * @param[in] massMatrix diagonal matrix M (stored as vector)
	 * @note Note that M^-1 * A is different from A * M^-1, even though M is diagonal
	 *       (M^-1*A: columns are multiplied by the same diag element)
	 */
	void divideByDiagonalMassMatrix(distributed_sparse_matrix& matrix,
			const distributed_vector& massMatrix);
	/**
	 * @short assemble and distribute internal face
	 * @param cell the cell to which the face belongs
	 * @param faceNumber the local number of the face (0,1,2 or 3)
	 * @param neighborCell the cell on the other side of the face
	 * @param neighborFaceNumber the local number of the face (0,1,2 or 3), as seen from the neigbor cell
	 * @param feFaceValues is passed in order to avoid allocating memory for each call of the function
	 * @param feSubfaceValues is passed in order to avoid allocating memory for each call of the function
	 * @param feNeighborFaceValues is passed in order to avoid allocating memory for each call of the function
	 */
	void assembleAndDistributeInternalFace(size_t direction,
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			size_t faceNumber,
			typename dealii::DoFHandler<dim>::cell_iterator& neighborCell,
			size_t neighborFaceNumber, dealii::FEFaceValues<dim>& feFaceValues,
			dealii::FESubfaceValues<dim>& feSubfaceValues,
			dealii::FEFaceValues<dim>& feNeighborFaceValues);

	/**
	 * @short map degrees of freedom to quadrature node indices on a cell
	 * @note called by the constructor to initialize m_dof_to_q_index
	 */
	std::map<size_t, size_t>  map_celldofs_to_q_index() const;

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
	 * @param[in] boltzmannModel the DQ model
	 */
	SEDGMinLee(shared_ptr<dealii::Triangulation<dim> > triangulation,
			shared_ptr<BoundaryCollection<dim> > boundaries,
			size_t orderOfFiniteElement,
			shared_ptr<BoltzmannModel> boltzmannModel, bool useCentralFlux = false);

	/// destructor
	virtual ~SEDGMinLee() {
		m_doFHandler->clear();
	}
	;

	/// function to (re-)assemble linear system
	virtual void reassemble();

	/// make streaming step
	virtual void stream();

	/// get global system matrix
	virtual const vector<distributed_sparse_matrix>& getSystemMatrix() const {
		return m_systemMatrix;
	}

	virtual void mapDoFsToSupportPoints(vector<dealii::Point<dim> >& supportPoints) const{
		dealii::DoFTools::map_dofs_to_support_points(m_mapping,
				*m_doFHandler, supportPoints);
	}


	virtual const shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const {
		return m_doFHandler;
	}

	const dealii::SparsityPattern& getSparsityPattern() const {
		return m_sparsityPattern;
	}

	const dealii::MappingQ1<dim>& getMapping() const {
		return m_mapping;
	}

	const std::map<size_t, size_t>& getCelldofToQIndex() const {
		return m_celldof_to_q_index;
	}

	const vector<std::map<size_t, size_t> >& getFacedofToQIndex() const {
		return m_facedof_to_q_index;
	}

	const shared_ptr<dealii::QGaussLobatto<dim - 1> >& getFaceQuadrature() const {
		return m_faceQuadrature;
	}

	const shared_ptr<dealii::FE_DGQArbitraryNodes<dim> >& getFe() const {
		return m_fe;
	}

	const shared_ptr<dealii::QGaussLobatto<dim> >& getQuadrature() const {
		return m_quadrature;
	}

	size_t getOrderOfFiniteElement() const {
		return m_orderOfFiniteElement;
	}
};

} /* namespace natrium */


#endif /* SEDGMINLEE_H_ */
