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
#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"

namespace natrium {


/* forward declaration */
class Stencil;

/**
 * @short Exception class for SEDGMinLee
 */
class AdvectionSolverException: public NATriuMException {
private:
	std::string message;
public:
	AdvectionSolverException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	AdvectionSolverException(const string& msg) :
			NATriuMException(msg), message(msg) {
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
 *         \partial_t f_{\alpha} + e_{\alpha} \partial_x f_{\alpha} = 0,\quad \forall {\alpha} = 1,\dots,Q-1
 *         \f]
 *         where \f$ f_{\alpha}(x,t) \f$ are the particle distribution functions, and \f$ e_{\alpha} \f$ are the particle
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
 *         \left( \partial_t f_{\alpha} + \partial_x (e_{\alpha} f_{\alpha}), \Phi \right)_{\Omega_e}
 *         = \left(n \left[ e_i f_{\alpha} - F^{\ast}_{\alpha}(f) \right], \Phi \right)_{\partial \Omega_e}.
 *         \f]
 *         In this formulation \f$ F^{\ast}_{i}(f) \f$ denotes the numerical fluxes. They can be be calculated
 *         as central fluxes or Lax-Friedrichs fluxes. Lax-Friedrichs is in general more accurate for the advection equation.
 *         For detailed information on the fluxes, see the cited paper.
 *         For spatial integration a Gauss-Lobatto quadrature is used, which has the advantage that the resulting mass matrix
 *         M_{\alpha} = (\psi_j, \psi_k)_{\Omega_e} is diagonal. This circumvents the solution of a linear equation system.
 *         Each advection equation leads to a ODE
 *         \f[ \partial_t f_{\alpha} = M_{\alpha}^{-1}(- e_{\alpha x} D_{{\alpha}x} - e_{{\alpha}y} D_{{\alpha}y} + R_{\alpha}) f_{\alpha} + B_i f_{{\alpha}^{\ast}} + b_{\alpha}.\f]
 *         Altogether, for the example of the D2Q9, the system becomes
 *         \f[ \partial_t f_{1,\dots,Q} =
 *         \left( \matrix{
 *         L_1 	& 	0 	& 	B_1 & 	0	& 	0 	& 	0	&	0	&	0 \cr
 *         0	&	L_2	&	0	&	B_2	&	0	&	0	&	0	&	0 \cr
 *         B_3	&	0	&	L_3	&	0	&	0	&	0	&	0	&	0 \cr
 *         0	&	B_4	&	0	&	L_4	&	0	&	0	&	0	&	0 \cr
 *         0 	& 	0	&	0	&	0	&	L_5 & 	0 	& 	B_5 & 	0 \cr
 *         0 	& 	0	&	0	&	0	&	0	&	L_6 & 	0 	& 	B_6\cr
 *         0 	& 	0	&	0	&	0	&	B_7 & 	0 	& 	L_7 & 	0 \cr
 *         0 	& 	0	&	0	&	0	&	0	&	B_8 & 	0 	& 	L_8
 *         } \right)
 *		   f_{1,\dots,Q}
 *		   +
 *		   \left( \matrix{
 *		   b_1  \cr b_2 \cr b_3 \cr b_4 \cr b_5 \cr b_6 \cr b_7 \cr b_8
 *		   }\right),
 *         \f]
 *         where \f$ L_{\alpha} = M_{\alpha}^{-1}(- e_{{\alpha}x} D_{{\alpha}x} - e_{{\alpha}y} D_{{\alpha}y} + R_{\alpha}) \f$.
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class SEDGMinLee: public AdvectionOperator<dim> {

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

#ifdef WITH_TRILINOS
	// locally owned degrees of freedom (for MPI parallelization)
	dealii::IndexSet m_locallyOwnedDofs;
	// locally relevant degrees of freedom (i.e. ghost layer cells)
	dealii::IndexSet m_locallyRelevantDofs;
#endif

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
			dealii::FEFaceValues<dim>& feNeighborFaceValues, const vector<double>& inverseLocalMassMatrix);

	/**
	 * @short calculate system diagonal block matrix  (Dx*eix + Dy*eiy)
	 */
	void calculateAndDistributeLocalStiffnessMatrix(size_t alpha,
			const vector<dealii::FullMatrix<double> > &derivativeMatrices,
			dealii::FullMatrix<double> &systemMatrix, const vector<double>& inverseLocalMassMatrix,
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
			dealii::FEFaceValues<dim>& feNeighborFaceValues, const vector<double>& inverseLocalMassMatrix);

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
	 */
	SEDGMinLee(boost::shared_ptr<Mesh<dim> > triangulation,
			boost::shared_ptr<BoundaryCollection<dim> > boundaries,
			size_t orderOfFiniteElement,
			boost::shared_ptr<Stencil> stencil,  bool useCentralFlux =
					false);

	/// destructor
	virtual ~SEDGMinLee() {
		m_doFHandler->clear();
	}
	;

	/// function to (re-)assemble linear system
	virtual void reassemble();

	virtual  void setupDoFs();

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

#ifndef WITH_TRILINOS
	/**
	 * @short get sparsity pattern
	 * @note not available for trilinos, as trilinos has an internal format for sparsity patterns (Epetra_CrsGraph)
	 */
	const dealii::SparsityPattern& getSparsityPattern(size_t i) const {
		return m_systemMatrix.block(i,i).get_sparsity_pattern();
	}
#endif

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

	virtual size_t getNumberOfDoFsPerCell() const{
		return m_fe->dofs_per_cell;
	}

	virtual const boost::shared_ptr<dealii::QGaussLobatto<dim> >& getQuadrature() const {
		return m_quadrature;
	}

	virtual size_t getOrderOfFiniteElement() const {
		return m_orderOfFiniteElement;
	}

	virtual size_t getNumberOfDoFs() const {
		return getSystemMatrix().block(0,0).n();
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

};

} /* namespace natrium */


#endif /* SEDGMINLEE_H_ */
