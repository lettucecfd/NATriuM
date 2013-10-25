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
#include "deal.II/fe/fe_q.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/lac/sparse_matrix.h"

#include "StreamingData.h"
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

	/// dealii::DoFHandler to distribute the degrees of freedom over the Triangulation
	dealii::DoFHandler<dim> m_doFHandler;

	/// dealii::FEValues object to get function values and derivatives at node points
	shared_ptr<dealii::FEValues<dim,dim> > m_feValues;

	/// Triangulation
	shared_ptr<dealii::Triangulation<dim> > m_tria;

	/// Finite Element function on one cell
	dealii::FE_Q<dim> m_fe;

	/// Sparsity Pattern of the sparse matrix
	dealii::SparsityPattern m_sparsityPattern;

	/// sparse system mass matrix (diagonal)
	// TODO this can be done more efficiently by storing the
	// (diagonal) mass matrix in a numeric_vector (and leaving
	// out the calculations for 0-entries)
	distributed_sparse_matrix m_massMatrix;

	/// System matrices
	vector<distributed_sparse_matrix> m_D;

	/// System
	numeric_vector m_systemRhs;



public:

	/// constructor
	/**
	 * @short Constructor
	 * @param[in] triangulation The global mesh.
	 */
	DataMinLee2011(
			shared_ptr<dealii::Triangulation<dim> > triangulation, size_t orderOfFiniteElement);

	/// destructor
	virtual ~DataMinLee2011(){};

	/// function to (re-)assemble linear system
	virtual void reassemble();
};

} /* namespace natrium */

#endif /* DATAMINLEE2011_H_ */
