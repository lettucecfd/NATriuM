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
#include "deal.II/lac/block_sparsity_pattern.h"
#include "deal.II/base/quadrature_lib.h"

#include "AdvectionOperator.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../boltzmannmodels/BoltzmannModel.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short Exception class for AdvectionSolver
 */
class AdvectionSolverException: public std::exception {
private:
	std::string message;
public:
	AdvectionSolverException(const char *msg) :
			message(msg) {
	}
	~AdvectionSolverException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/** @short This class solves the linear advection equations by a scheme which is used, e.g.,
 * 		   by Min and Lee (2011): A spectral-element discontinuous Galerkin lattice Boltzmann
 * 		   method for nearly incompressible flows, JCP 230 pp. 245-259.
 *         The advection equations used in the Lattice Boltzmann on unstructured grids are
 *         \f[
 *         \partial_t f_i + e_i \partial_x f_i = 0,\quad \forall i = 1,\dots,Q-1
 *         \f]
 *         where \f$ f_i(x,t) \f$ are the particle distribution functions, and \f$ e_i \f$ are the particle
 *         velocities. The discontinuous Galerkin (DG) method turns these PDEs into a large system of ODEs which can then be solved
 *         by a time integration scheme. Whereas this class implements the SEDG spatial discretization, the
 *         time integration is done by a subclass of TimeIntegrator, e.g. RungeKutta5LowStorage.
 *         In other Finite Element schemes, degrees of freedom can belong to different elements
 *         (e.g. at corners of elements). In contrast, DG methods have the degrees of freedom belonging to a
 *         single element, which can lead to discontinuities at the element faces. To connect neighbor cells,
 *         the integral over the boundary of each cell incorporates the solution on neighbor cells. These
 *         contributions are called numerical fluxes.
 *         The DG scheme uses the weak formulation of the above equations on quadrilateral elements $\Omega_e$:
 *         \f[
 *         \left( \partial_t f_i + \partial_x (e_i f_i), \Phi \right)_{\Omega_e}
 *         = \left(n \left[ e_i f_i - F^{\ast}_{i}(f) \right], \Phi \right)_{\partial \Omega_e}.
 *         \f]
 *         In this formulation \f$ F^{\ast}_{i}(f) \f$ denotes the numerical fluxes. They can be be calculated
 *         as central fluxes or Lax-Friedrichs fluxes. Lax-Friedrichs is in general more accurate for the advection equation.
 *         For detailed information on the fluxes, see the cited paper.
 *         For spatial integration a Gauss-Lobatto quadrature is used, which has the advantage that the resulting mass matrix
 *         M_i = (\psi_j, \psi_k)_{\Omega_e} is diagonal. This circumvents the solution of a linear equation system.
 *         Each advection equation leads to a ODE
 *         \f[ \partial_t f_i = M_i^{-1}(- e_{ix} D_{ix} - e_{iy} D_{iy} + R_i) f_i + B_i f_{i^{\ast}} + b_i.\f]
 *         Altogether, for the example of the D2Q9, the system becomes
 *         \f[
 *
 *         \f]
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
	dealii::BlockSparsityPattern m_sparsityPattern;

	/// Mapping from real space to unit cell
	const dealii::MappingQ1<dim> m_mapping;

	/// global mass matrix M
	distributed_vector m_massMatrix;

	/// System matrix L = M^(-1)*(-D+R)
	distributed_sparse_block_matrix m_systemMatrix;

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
	void divideByDiagonalMassMatrix(dealii::SparseMatrix<double>& matrix,
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

	/**
	 * calculate a 'magic number' which helps to test whether two triangulations are equal
	 */
	double calcMagicNumber() const;

	/** @short save matrices to files
	 *  @param[in] directory directory to save the matrix files to
	 */
	void saveMatricesToFiles(const string& directory) const;

	/** @short load matrices from files
	 *  @param[in] directory directory to load the matrix files from
	 */
	void loadMatricesFromFiles(const string& directory);

	/**
	 * @short A function to write the status of the solver, mesh, etc at checkpoint time in order to prevent corrupt restarts of a simulation
	 * @param[in] directory the output directory of a simulation
	 */
	void writeStatus(const string& directory) const;

	/**
	 * @short A function to write the status of the solver, mesh, etc at checkpoint time in order to prevent corrupt restarts of a simulation
	 * @param[in] directory the output directory of a simulation
	 * @param[out] if not OK, the message indicates what is not OK
	 */
	bool isStatusOK(const string& directory, string& message) const;

	/** @short load matrices and status from files
	 *  @param[in] directory directory to load the matrix files from
	 *  @throws AdvectionSolverException
	 */
	void loadCheckpoint(const string& directory);

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
			shared_ptr<BoltzmannModel> boltzmannModel, string inputDirectory = "", bool useCentralFlux =
					false);

	/// destructor
	virtual ~SEDGMinLee() {
		m_doFHandler->clear();
	}
	;

	/// function to (re-)assemble linear system
	virtual void reassemble();

	/// make streaming step
	virtual void stream();

	/** @short save matrices and status to files
	 *  @param[in] directory directory to save the matrix files to
	 *  @throws AdvectionSolverException
	 */
	virtual void saveCheckpoint(const string& directory) const;

	/// get global system matrix
	virtual const distributed_sparse_block_matrix& getSystemMatrix() const {
		return m_systemMatrix;
	}

	virtual void mapDoFsToSupportPoints(
			vector<dealii::Point<dim> >& supportPoints) const {
		assert(supportPoints.size() == this->getNumberOfDoFs());
		dealii::DoFTools::map_dofs_to_support_points(m_mapping, *m_doFHandler,
				supportPoints);
	}

	virtual const shared_ptr<dealii::DoFHandler<dim> >& getDoFHandler() const {
		return m_doFHandler;
	}

	const dealii::BlockSparsityPattern& getBlockSparsityPattern() const {
		return m_sparsityPattern;
	}

	const dealii::SparsityPattern& getSparsityPattern(size_t i) const {
		return m_systemMatrix.block(i,i).get_sparsity_pattern();
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

	virtual size_t getNumberOfDoFs() const {
		return getSystemMatrix().block(0,0).n();
	}
};

} /* namespace natrium */


#endif /* SEDGMINLEE_H_ */
