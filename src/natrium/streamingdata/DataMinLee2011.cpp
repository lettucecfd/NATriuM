/**
 * @file DataMinLee2011.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DataMinLee2011.h"

#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
DataMinLee2011<dim>::DataMinLee2011(
		shared_ptr<Triangulation<dim> > triangulation,
		shared_ptr<BoundaryCollection<dim> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel) :
		m_tria(triangulation), m_boundaries(boundaries), m_boltzmannModel(
				boltzmannModel) {
	// make dof handler
	/*	m_quadrature = make_shared<QGaussLobatto<dim> >(orderOfFiniteElement);
	 m_fe = make_shared<FE_DGQArbitraryNodes<dim> >(*m_quadrature);
	 m_doFHandler = make_shared<DoFHandler<dim> >(*triangulation);
	 */
	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);
	/*
	 // make sparsity pattern
	 updateSparsityPattern();
	 updateSystemMatrixAccordingToSparsityPattern();
	 */
	//make sparse matrix
	/*CompressedSparsityPattern cSparse(m_doFHandler->n_dofs());

	 //reorder degrees of freedom
	 DoFRenumbering::Cuthill_McKee(*m_doFHandler);
	 DoFTools::make_sparsity_pattern(*m_doFHandler, cSparse);
	 m_sparsityPattern.copy_from(cSparse);

	 //reinitialize matrices
	 m_massMatrix.reinit(m_sparsityPattern);
	 for (size_t i = 0; i < dim; i++) {
	 distributed_sparse_matrix D_i;
	 m_derivativeMatrix.push_back(D_i);
	 m_derivativeMatrix.at(i).reinit(m_sparsityPattern);
	 }

	 for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
	 distributed_sparse_matrix R_i;
	 m_faceMatrix.push_back(R_i);
	 m_faceMatrix.at(i).reinit(m_sparsityPattern);
	 }

	 for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
	 distributed_sparse_matrix L;
	 m_systemMatrix.push_back(L);
	 }

	 // fe values
	 m_feValues = make_shared<FEValues<dim, dim> >(m_fe, m_quadrature,
	 update_values | update_gradients | update_JxW_values);

	 // assemble system
	 reassemble();
	 */
} /* DataMinLee2011<dim>::DataMinLee2011 */
/// The template parameter must be made explicit in order for the code to compile
template DataMinLee2011<2>::DataMinLee2011(
		shared_ptr<Triangulation<2> > triangulation,
		shared_ptr<BoundaryCollection<2> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel);
template DataMinLee2011<3>::DataMinLee2011(
		shared_ptr<Triangulation<3> > triangulation,
		shared_ptr<BoundaryCollection<3> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel);

template<size_t dim>
inline void DataMinLee2011<dim>::updateSparsityPattern() {

	//make sparse matrix
	CompressedSparsityPattern cSparse(m_doFHandler->n_dofs());

	//reorder degrees of freedom
	DoFRenumbering::Cuthill_McKee(*m_doFHandler);
	DoFTools::make_flux_sparsity_pattern(*m_doFHandler, cSparse);

	size_t dofs_per_cell = m_doFHandler->get_fe().dofs_per_cell;

	//add periodic neighbors
	const vector<shared_ptr<PeriodicBoundary1D> > periodicBoundaries =
			m_boundaries->getPeriodicBoundaries();
	for (size_t i = 0; i < periodicBoundaries.size(); i++) {
		// map cells to each other
		// TODO only update sparsity pattern for changed cells
		periodicBoundaries.at(i)->createCellMap(*m_doFHandler);
		const std::map<dealii::DoFHandler<2>::active_cell_iterator,
				std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t> > cellMap =
				periodicBoundaries.at(i)->getCellMap();
		std::map<dealii::DoFHandler<2>::active_cell_iterator,
				std::pair<dealii::DoFHandler<2>::active_cell_iterator, size_t> >::const_iterator element =
				cellMap.begin();
		// for each cells belonging to the periodic boundary
		for (; element != cellMap.end(); element++) {
			vector<dealii::types::global_dof_index> doFIndicesAtCell1(dofs_per_cell);
			vector<dealii::types::global_dof_index> doFIndicesAtCell2(dofs_per_cell);
			element->first->get_dof_indices(doFIndicesAtCell1);
			element->first->get_dof_indices(doFIndicesAtCell2);
			// couple all dofs at boundary 1 with dofs at boundary 2
			// TODO only couple the ones which are nonzero at the face (are there any???)
			for (size_t j = 0; j < dofs_per_cell; j++) {
				for (size_t k = 0; k < dofs_per_cell; k++) {
					cSparse.add(j,k);
				}
			}
		}
	}

	m_sparsityPattern.copy_from(cSparse);

	//reinitialize matrices
	for (size_t i = 0; i < m_systemMatrix.size(); i++) {
		m_systemMatrix.at(i).reinit(m_sparsityPattern);
	}
} /* updateSparsityPattern */
// The template parameter has to be made expicit in order for the code to compile
template void DataMinLee2011<2>::updateSparsityPattern();
//template void DataMinLee2011<3>::updateSparsityPattern();

template<size_t dim>
void DataMinLee2011<dim>::stream() {
}
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::stream();
template void DataMinLee2011<3>::stream();

template<size_t dim>
void DataMinLee2011<dim>::reassemble() {
	// TODO: if Triangulation changed: reinit dof-handler and sparsity pattern in some way
	/*
	 const unsigned int dofs_per_cell = m_fe->dofs_per_cell;
	 const unsigned int n_q_points = m_quadrature->size();

	 // Initialize matrices
	 dealii::FullMatrix<double> localMassMatrix;
	 vector<dealii::FullMatrix<double> > localDerivativeMatrix;
	 for (size_t i = 0; i < dim; i++) {
	 dealii::FullMatrix<double> D_i(dofs_per_cell, dofs_per_cell);
	 localDerivativeMatrix.push_back(D_i);
	 }
	 dealii::FullMatrix<double> localFaceMatrix(dofs_per_cell);
	 dealii::FullMatrix<double> localSystemMatrix(dofs_per_cell);
	 std::vector<types::global_dof_index> localDofIndices(dofs_per_cell);

	 ///////////////
	 // MAIN LOOP //
	 ///////////////
	 DoFHandler<2>::active_cell_iterator cell = m_doFHandler->begin_active(),
	 endc = m_doFHandler->end();
	 for (; cell != endc; ++cell) {
	 // initialize
	 m_feValues->reinit(cell);
	 localMassMatrix = 0;
	 localDerivativeMatrix = 0;
	 localFaceMatrix = 0;
	 localSystemMatrix = 0;

	 // assemble mass matrix
	 for (size_t i = 0; i < dofs_per_cell; ++i)
	 for (size_t j = 0; j < dofs_per_cell; ++j)
	 for (size_t q_point = 0; q_point < n_q_points; ++q_point)
	 localMassMatrix(i, j) += (m_feValues->shape_value(i, q_point)
	 * m_feValues->value(j, q_point)
	 * m_feValues->JxW(q_point));

	 // assemble derivative matrices
	 for (size_t i = 0; i < dim; i++) {
	 for (size_t i = 0; i < dofs_per_cell; ++i)
	 for (size_t j = 0; j < dofs_per_cell; ++j)
	 for (size_t q_point = 0; q_point < n_q_points; ++q_point)
	 localDerivativeMatrix(i, j) +=
	 (m_feValues->shape_gradient(i, q_point)
	 * m_feValues->shape_gradient(j, q_point)
	 * m_feValues->JxW(q_point));
	 }

	 for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
	 // collect all neighbors
	 // TODO What about hanging nodes with DG methods.


	 // assemble face matrix
	 for (size_t i = 0; i < dofs_per_cell; ++i)
	 for (size_t j = 0; j < dofs_per_cell; ++j)
	 for (size_t q_point = 0; q_point < n_q_points; ++q_point)
	 localFaceMatrix(i, j) +=
	 (m_feValues->normal_vector(i, q_point)
	 * m_feValues->shape_gradient(j, q_point)
	 * m_feValues->JxW(q_point));
	 }
	 }
	 */
}
/// The template parameter must be made explicit in order for the code to compile
template void DataMinLee2011<2>::reassemble();
template void DataMinLee2011<3>::reassemble();

} /* namespace natrium */
