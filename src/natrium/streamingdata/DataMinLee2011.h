/**
 * @file DataMinLee2011.h
 * @short Global data which is used by Min and Lee (2011): A spectral-elemennt discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef DATAMINLEE2011_H_
#define DATAMINLEE2011_H_

#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/base/quadrature_lib.h"


#include "StreamingData.h"
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
template <size_t dim> class DataMinLee2011: public  StreamingData<dim>{

private:

	/// Triangulation
	shared_ptr<dealii::Triangulation<dim> > m_tria;

	/// Boundary Description
	shared_ptr<BoundaryCollection<dim> > m_boundaries;

	/// integration on gauss lobatto nodes
	shared_ptr<dealii::QGaussLobatto<dim> > m_quadrature;

	/// Finite Element function on one cell
	shared_ptr<dealii::FE_DGQArbitraryNodes<dim> > m_fe;

	/// dealii::DoFHandler to distribute the degrees of freedom over the Triangulation
	shared_ptr<dealii::DoFHandler<dim> > m_doFHandler;

	/// dealii::FEValues object to get function values and derivatives at node points
	shared_ptr<dealii::FEValues<dim,dim> > m_feValues;

	/// Sparsity Pattern of the sparse matrix
	dealii::SparsityPattern m_sparsityPattern;

	/// System matrix L = M^(-1)*(-D+R)
	vector<distributed_sparse_matrix> m_systemMatrix;

	/// the DQ model (e.g. D2Q9)
	shared_ptr<BoltzmannModel> m_boltzmannModel;

	/// number of dofs per cell
	size_t m_dofs_per_cell;

	// number of integration points
	size_t m_n_quadrature_points;

	/**
	 * @short update the sparsity pattern of the system matrix
	 */
	void updateSparsityPattern();

	/**
	 * @short update the global system matrix
	 */
	void updateSystemMatrixAccordingToSparsityPattern();

	/**
	 * @short assemble local mass matrix M
	 * @param[out] massMatrix The mass matrix <phi_i, phi_j>
	 *             For SEDG methods the mass Matrix is diagonal.
	 *
	 */
	// TODO SEDG implemenation with fully diagonal mass matrix
	void assembleLocalMassMatrix(dealii::FullMatrix<double> &massMatrix) const;

	/**
	 * @short assemble the i-th local derivative matrix
	 * @param[in] i < dim; 0 for Dx, 1 for Dy, 2 for Dz (3D)
	 * @param[out] derivativeMatrix The i-th derivative matrix <D_i phi_j, phi_k>
	 */
	void assembleLocalDerivativeMatrix(size_t i, dealii::FullMatrix<double> &derivativeMatrix) const;

	/**
	 * @short assemble local face matrix;
	 * @param[in] i < Q; this matrix is dependent on the flow direction
	 * @param[out] faceMatrix The integral over all faces, incorporating boundary conditions
	 */
	void assembleLocalFaceMatrix(size_t i, dealii::FullMatrix<double> &faceMatrix) const;

	/**
	 * @short calculate system matrix L
	 */
	void calculateLocalSystemMatrix();


public:

	/// constructor
	/**
	 * @short Constructor
	 * @param[in] triangulation The global mesh.
	 * @param[in] orderOfFiniteElement The number of nodes element and dimension
	 * @param[in] boltzmannModel the DQ model
	 */
	DataMinLee2011(
			shared_ptr<dealii::Triangulation<dim> > triangulation, shared_ptr<BoundaryCollection<dim> > boundaries, size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel);

	/// destructor
	virtual ~DataMinLee2011(){};

	/// function to (re-)assemble linear system
	virtual void reassemble();

	/// make streaming step
	virtual void stream();
};

} /* namespace natrium */


#endif /* DATAMINLEE2011_H_ */
