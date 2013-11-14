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
	m_quadrature = make_shared<QGaussLobatto<dim> >(orderOfFiniteElement);
	m_fe = make_shared<FE_DGQArbitraryNodes<dim> >(QGaussLobatto<1>(orderOfFiniteElement));
	m_doFHandler = make_shared<DoFHandler<dim> >(*triangulation);
	m_dofs_per_cell = m_fe->dofs_per_cell;
	m_n_quadrature_points = m_quadrature->size();

	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);
	updateSparsityPattern();

	 // fe values
	 m_feValues = make_shared<FEValues<dim, dim> >(*m_fe, *m_quadrature,
	 update_values | update_gradients | update_JxW_values);

	 // assemble system
	 reassemble();
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
	const vector<shared_ptr<PeriodicBoundary<dim> > > periodicBoundaries =
			m_boundaries->getPeriodicBoundaries();
	for (size_t i = 0; i < periodicBoundaries.size(); i++) {
		// map cells to each other
		// TODO only update sparsity pattern for changed cells
		periodicBoundaries.at(i)->createCellMap(*m_doFHandler);
		const std::map<typename dealii::DoFHandler<dim>::active_cell_iterator,
				std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, size_t> > cellMap =
				periodicBoundaries.at(i)->getCellMap();
		typename std::map<typename dealii::DoFHandler<dim>::active_cell_iterator,
				std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, size_t> >::const_iterator element =
				cellMap.begin();
		// for each cells belonging to the periodic boundary
		for (; element != cellMap.end(); element++) {
			vector<dealii::types::global_dof_index> doFIndicesAtCell1(
					dofs_per_cell);
			vector<dealii::types::global_dof_index> doFIndicesAtCell2(
					dofs_per_cell);
			element->first->get_dof_indices(doFIndicesAtCell1);
			element->first->get_dof_indices(doFIndicesAtCell2);
			// couple all dofs at boundary 1 with dofs at boundary 2
			// TODO only couple the ones which are nonzero at the face (are there any???)
			// TODO remove the INVARIANT "discretization at boundary 1 = discretization at boundary 2"
			//      e.g. by mapping, allowing more than one periodic neighbor, ...
			for (size_t j = 0; j < dofs_per_cell; j++) {
				for (size_t k = 0; k < dofs_per_cell; k++) {
					cSparse.add(j, k);
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
//TODO generalize to 3D
//template void DataMinLee2011<3>::updateSparsityPattern();



template<size_t dim>
inline void DataMinLee2011<dim>::assembleLocalMassMatrix(
		dealii::FullMatrix<double>& massMatrix) const {
	massMatrix = 0;
	for (size_t i = 0; i < m_dofs_per_cell; ++i)
		for (size_t j = 0; j < m_dofs_per_cell; ++j)
			for (size_t q_point = 0; q_point < m_n_quadrature_points; ++q_point)
				massMatrix(i, j) += (m_feValues->shape_value(i,
						q_point) * m_feValues->shape_value(j, q_point)
						* m_feValues->JxW(q_point));
} /*assembleLocalMassMatrix*/
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::assembleLocalMassMatrix(
		dealii::FullMatrix<double>& massMatrix) const;
template void DataMinLee2011<3>::assembleLocalMassMatrix(
		dealii::FullMatrix<double>& massMatrix) const;


template<size_t dim>
inline void DataMinLee2011<dim>::assembleLocalDerivativeMatrix(size_t i,
		dealii::FullMatrix<double>& derivativeMatrix) const {
	derivativeMatrix = 0;
	for (size_t i = 0; i < m_dofs_per_cell; ++i)
		for (size_t j = 0; j < m_dofs_per_cell; ++j)
			for (size_t q_point = 0; q_point < m_n_quadrature_points; ++q_point)
				derivativeMatrix(i, j) +=
						(m_feValues->shape_grad(i, q_point)
								* m_feValues->shape_grad(j, q_point)
								* m_feValues->JxW(q_point));
} /* assembleLocalDerivativeMatrix */
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::assembleLocalDerivativeMatrix(size_t i,
		dealii::FullMatrix<double>& derivativeMatrix) const;
template void DataMinLee2011<3>::assembleLocalDerivativeMatrix(size_t i,
		dealii::FullMatrix<double>& derivativeMatrix) const;


template<size_t dim>
void DataMinLee2011<dim>::stream() {
}
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::stream();
template void DataMinLee2011<3>::stream();

template<size_t dim>
void DataMinLee2011<dim>::reassemble() {
	// TODO: if Triangulation changed: reinit dof-handler and sparsity pattern in some way

	// Initialize matrices
	dealii::FullMatrix<double> localMassMatrix;
	vector<dealii::FullMatrix<double> > localDerivativeMatrices;
	for (size_t i = 0; i < dim; i++) {
		dealii::FullMatrix<double> D_i(m_dofs_per_cell, m_dofs_per_cell);
		localDerivativeMatrices.push_back(D_i);
	}
	dealii::FullMatrix<double> localFaceMatrix(m_dofs_per_cell);
	dealii::FullMatrix<double> localSystemMatrix(m_dofs_per_cell);
	std::vector<types::global_dof_index> localDofIndices(m_dofs_per_cell);

	///////////////
	// MAIN LOOP //
	///////////////
	typename DoFHandler<dim>::active_cell_iterator cell = m_doFHandler->begin_active(),
			endc = m_doFHandler->end();
	for (; cell != endc; ++cell) {
		// calculate the fe values for the cell
		m_feValues->reinit(cell);

		// assemble local matrices
		assembleLocalMassMatrix(localMassMatrix);
		for (size_t i = 0; i < dim; i++){
			assembleLocalDerivativeMatrix(i, localDerivativeMatrices.at(i));
		}
		for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
			assembleLocalFaceMatrix(i, localFaceMatrix);

			// calculate local system matrix L
			//calculateLocalSystemMatrix(localMassMatrix, localDerivativeMatrices, localFaceMatrix, localSystemMatrix);

			//

			// collect all neighbors
			// TODO What about hanging nodes with DG methods.

			// assemble face matrix
			/*
			 *
			 *for (size_t i = 0; i < dofs_per_cell; ++i)
				for (size_t j = 0; j < dofs_per_cell; ++j)
					for (size_t q_point = 0; q_point < n_q_points; ++q_point)
						localFaceMatrix(i, j) += (m_feValues->normal_vector(i,
								q_point)
								* m_feValues->shape_gradient(j, q_point)
								* m_feValues->JxW(q_point));
		*/}
	}

}
/// The template parameter must be made explicit in order for the code to compile
template void DataMinLee2011<2>::reassemble();
template void DataMinLee2011<3>::reassemble();

} /* namespace natrium */
