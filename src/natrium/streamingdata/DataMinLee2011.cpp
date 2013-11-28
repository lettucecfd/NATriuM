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
#include "deal.II/fe/fe_update_flags.h"

#include "../problemdescription/PeriodicBoundary.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
DataMinLee2011<dim>::DataMinLee2011(
		shared_ptr<Triangulation<dim> > triangulation,
		shared_ptr<BoundaryCollection<dim> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel) :
		m_tria(triangulation), m_boundaries(boundaries), m_mapping(), m_boltzmannModel(
				boltzmannModel) {
	// make dof handler
	m_quadrature = make_shared<QGaussLobatto<dim> >(orderOfFiniteElement);
	m_faceQuadrature = make_shared<QGaussLobatto<dim - 1> >(
			orderOfFiniteElement);
	m_fe = make_shared<FE_DGQArbitraryNodes<dim> >(
			QGaussLobatto<1>(orderOfFiniteElement));
	m_doFHandler = make_shared<DoFHandler<dim> >(*triangulation);

	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);
	updateSparsityPattern();

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
void DataMinLee2011<dim>::updateSparsityPattern() {

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
				std::pair<typename dealii::DoFHandler<dim>::cell_iterator,
						size_t> > cellMap =
				periodicBoundaries.at(i)->getCellMap();
		typename std::map<
				typename dealii::DoFHandler<dim>::active_cell_iterator,
				std::pair<typename dealii::DoFHandler<dim>::cell_iterator,
						size_t> >::const_iterator element = cellMap.begin();
		// for each cells belonging to the periodic boundary
		for (; element != cellMap.end(); element++) {
			vector<dealii::types::global_dof_index> doFIndicesAtCell1(
					dofs_per_cell);
			vector<dealii::types::global_dof_index> doFIndicesAtCell2(
					dofs_per_cell);
			element->first->get_dof_indices(doFIndicesAtCell1);
			element->second.first->get_dof_indices(doFIndicesAtCell2);
			// couple all dofs at boundary 1 with dofs at boundary 2
			// TODO only couple the ones which are nonzero at the face (are there any???)
			// TODO remove the INVARIANT "discretization at boundary 1 = discretization at boundary 2"
			//      e.g. by mapping, allowing more than one periodic neighbor, ...
			for (size_t j = 0; j < dofs_per_cell; j++) {
				for (size_t k = 0; k < dofs_per_cell; k++) {
					cSparse.add(doFIndicesAtCell1.at(j),
							doFIndicesAtCell2.at(k));
				}
			}
		}
	}

	m_sparsityPattern.copy_from(cSparse);

	//reinitialize matrices
	for (size_t i = 0; i < m_systemMatrix.size(); i++) {
		m_systemMatrix.at(i).reinit(m_sparsityPattern);
	}

	m_massMatrix.reinit(m_doFHandler->n_dofs());
	std::fill(m_massMatrix.begin(), m_massMatrix.end(), 0.0);

} /* updateSparsityPattern */
// The template parameter has to be made expicit in order for the code to compile
template void DataMinLee2011<2>::updateSparsityPattern();
//TODO generalize to 3D
//template void DataMinLee2011<3>::updateSparsityPattern();

template<size_t dim>
void DataMinLee2011<dim>::assembleLocalMassMatrix(
		const dealii::FEValues<dim>& feValues, size_t dofs_per_cell,
		size_t n_q_points, vector<double> &massMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs) {
	// initialize with zeros
	// TODO the fill operation can be cut out in the final implementation,
	// changing the += in the loop to =
	std::fill(massMatrix.begin(), massMatrix.end(), 0.0);
	// fill diagonal "matrix"
	for (size_t i = 0; i < dofs_per_cell; ++i)
		// TODO the inner for loop can be cut out in final implementation
		for (size_t j = 0; j < dofs_per_cell; ++j) {
			// TODO nonDiagonalElement is just for testing; cut out in final implementation
			double nonDiagonalElement = 0.0;
			// TODO the most inner loop can be cut out in the final implementation as the
			// shape values are 0 for all nodes except the base node
			for (size_t q_point = 0; q_point < n_q_points; ++q_point) {
				if (i != j) {
					nonDiagonalElement += feValues.shape_value(i, q_point)
							* feValues.shape_value(j, q_point)
							* feValues.JxW(q_point);
				} else {
					massMatrix.at(i) += feValues.shape_value(i, q_point)
							* feValues.shape_value(j, q_point)
							* feValues.JxW(q_point);
				}
			}
			assert(abs(nonDiagonalElement) < 1e-15);
		}
	// Assemble to global mass matrix
	for (size_t i = 0; i < dofs_per_cell; i++) {
		m_massMatrix(globalDoFs.at(i)) = massMatrix.at(i);
	}
} /*assembleLocalMassMatrix*/
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::assembleLocalMassMatrix(
		const dealii::FEValues<2>& feValues, size_t dofs_per_cell,
		size_t n_q_points, vector<double> &massMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs);
template void DataMinLee2011<3>::assembleLocalMassMatrix(
		const dealii::FEValues<3>& feValues, size_t dofs_per_cell,
		size_t n_q_points, vector<double> &massMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs);

template<size_t dim>
void DataMinLee2011<dim>::assembleLocalDerivativeMatrix(size_t i,
		const dealii::FEValues<dim>& feValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double> &derivativeMatrix) const {
	derivativeMatrix = 0;
	for (size_t i = 0; i < dofs_per_cell; ++i)
		for (size_t j = 0; j < dofs_per_cell; ++j)
			for (size_t q_point = 0; q_point < n_q_points; ++q_point)
				derivativeMatrix(i, j) += (feValues.shape_grad(i, q_point)
						* feValues.shape_grad(j, q_point)
						* feValues.JxW(q_point));
} /* assembleLocalDerivativeMatrix */
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::assembleLocalDerivativeMatrix(size_t i,
		const dealii::FEValues<2>& feValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double> &derivativeMatrix) const;
template void DataMinLee2011<3>::assembleLocalDerivativeMatrix(size_t i,
		const dealii::FEValues<3>& feValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double> &derivativeMatrix) const;

template<size_t dim>
void DataMinLee2011<dim>::assembleAndDistributeLocalFaceMatrices(size_t i,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		dealii::FEFaceValuesBase<dim>& feFaceValues,
		dealii::FEFaceValuesBase<dim>& feSubfaceValues,
		dealii::FEFaceValuesBase<dim>& feNeighborFaceValues,
		size_t dofs_per_cell, size_t n_q_points,
		dealii::FullMatrix<double>& faceMatrix) {
	bool assemblyDone = false;
	// loop over all faces
	for (size_t j = 0; j < dealii::GeometryInfo<2>::faces_per_cell; j++) {
		assemblyDone = false;

		//Faces at boundary
		if (cell->face(j)->at_boundary()) {
			// TODO At this point the implementation is not efficient: Loop over all Boundaries :/ pfui
			// Better: BoundaryCollection should have a map < <Boundary cell, face_id>, Boundary ID>

			// Apply periodic boundaries
			for (size_t k = 0; k < m_boundaries->numberOfPeriodicBoundaries();
					k++) {
				const shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
						m_boundaries->getPeriodicBoundaries().at(k);
				if (periodicBoundary->isFaceInBoundary(cell, j)) {
					typename dealii::DoFHandler<dim>::cell_iterator neighborCell;
					periodicBoundary->getOppositeCellAtPeriodicBoundary(cell,
							neighborCell);
					assembleAndDistributeInternalFace(
							cell, j, neighborCell,
							dealii::GeometryInfo<dim>::opposite_face[j],
							feFaceValues, feSubfaceValues,
							feNeighborFaceValues);
					assemblyDone = true;
					break;
				}
			} /* for all periodic boundaries */
			if (assemblyDone) {
				continue;
			}

			// Apply other boundaries
			// TODO Implement other boundary conditions
		} else {
			// Internal faces
			typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(
					j);
			assembleAndDistributeInternalFace(cell, j,
					neighbor, dealii::GeometryInfo<dim>::opposite_face[j],
					feFaceValues, feSubfaceValues, feNeighborFaceValues);
		} /* if (face at boundary) {} else {} */

	}

} /* assembleLocalFaceMatrix */
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::assembleAndDistributeLocalFaceMatrices(
		size_t i,
		typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		dealii::FEFaceValuesBase<2>& feFaceValues,
		dealii::FEFaceValuesBase<2>& feSubfaceValues,
		dealii::FEFaceValuesBase<2>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix);
template void DataMinLee2011<3>::assembleAndDistributeLocalFaceMatrices(
		size_t i,
		typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		dealii::FEFaceValuesBase<3>& feFaceValues,
		dealii::FEFaceValuesBase<3>& feSubfaceValues,
		dealii::FEFaceValuesBase<3>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix);

template<size_t dim>
void DataMinLee2011<dim>::divideByDiagonalMassMatrix(
		distributed_sparse_matrix& matrix,
		const distributed_vector& massMatrix) {
	size_t n = massMatrix.size();
	for (size_t i; i < n; i++) {
		for (size_t j; j < n; j++)
			// TODO Parallel loop
			if (m_sparsityPattern.exists(i, j))
				matrix.set(i, j, matrix(i, j) / massMatrix(j));
	}
}
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::divideByDiagonalMassMatrix(
		distributed_sparse_matrix& matrix,
		const distributed_vector& massMatrix);
template void DataMinLee2011<3>::divideByDiagonalMassMatrix(
		distributed_sparse_matrix& matrix,
		const distributed_vector& massMatrix);

template<> void DataMinLee2011<2>::calculateAndDistributeLocalStiffnessMatrix(
		size_t i, const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs,
		size_t dofsPerCell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	systemMatrix = derivativeMatrices.at(0);
	systemMatrix *= (-m_boltzmannModel->getDirection(i)[0]);
	systemMatrix.add(-m_boltzmannModel->getDirection(i)[1],
			derivativeMatrices.at(1));
// distribute to global system matrix
	for (unsigned int j = 0; i < dofsPerCell; j++)
		for (unsigned int k = 0; j < dofsPerCell; k++)
			m_systemMatrix.at(i).add(globalDoFs[j], globalDoFs[k],
					systemMatrix(j, k));
}
template<> void DataMinLee2011<3>::calculateAndDistributeLocalStiffnessMatrix(
		size_t i, const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs,
		size_t dofsPerCell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	systemMatrix = derivativeMatrices.at(0);
	systemMatrix *= (-m_boltzmannModel->getDirection(i)[0]);
	systemMatrix.add(-m_boltzmannModel->getDirection(i)[1],
			derivativeMatrices.at(1), -m_boltzmannModel->getDirection(i)[2],
			derivativeMatrices.at(2));
// distribute to global system matrix
	for (unsigned int j = 0; i < dofsPerCell; j++)
		for (unsigned int k = 0; j < dofsPerCell; k++)
			m_systemMatrix.at(i).add(globalDoFs[j], globalDoFs[k],
					systemMatrix(j, k));
}

template<size_t dim>
void DataMinLee2011<dim>::assembleAndDistributeInternalFace(
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<dim>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValuesBase<dim>& feFaceValues,
		dealii::FEFaceValuesBase<dim>& feSubfaceValues,
		dealii::FEFaceValuesBase<dim>& feNeighborFaceValues) {
// get the required FE Values for the local cell
	//feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<Point<dim> > &normals = feFaceValues.get_normal_vectors();

// get the required dofs of the neighbor cell
	typename DoFHandler<dim>::face_iterator neighborFace = neighborCell->face(
			neighborFaceNumber);

// WARNING: Hanging nodes are not implemented, yet
// TODO Implement local refinement
	/*if (not neighborFace->has_children()){
	 if () {

	 }
	 feNeighborFaceValues.reinit(neighborCell, neighborFaceNumber);
	 }*/
// TODO efficient multiplication with unit vector
	vector<Point<dim> > beta(feFaceValues.n_quadrature_points);
//beta_function.value_list (feFaceValues.get_quadrature_points(), beta);
	for (unsigned int point = 0; point < feFaceValues.n_quadrature_points;
			++point) {
		const double beta_n = beta[point] * normals[point];
		if (beta_n > 0) {
		}
	}

} /* assembleAndDistributeInternalFace */
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::assembleAndDistributeInternalFace(
		typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<2>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValuesBase<2>& feFaceValues,
		dealii::FEFaceValuesBase<2>& feSubfaceValues,
		dealii::FEFaceValuesBase<2>& feNeighborFaceValues);
template void DataMinLee2011<3>::assembleAndDistributeInternalFace(
		typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<3>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValuesBase<3>& feFaceValues,
		dealii::FEFaceValuesBase<3>& feSubfaceValues,
		dealii::FEFaceValuesBase<3>& feNeighborFaceValues);

template<size_t dim>
void DataMinLee2011<dim>::stream() {
}
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::stream();
template void DataMinLee2011<3>::stream();

template<size_t dim>
void DataMinLee2011<dim>::reassemble() {
// TODO: if Triangulation changed: reinit dof-handler and sparsity pattern in some way

/////////////////////////////////
// Initialize Finite Element ////
/////////////////////////////////
// Define update flags (which values have to be known at each cell, face, neighbor face)
	const dealii::UpdateFlags cellUpdateFlags = update_values | update_gradients
			| update_quadrature_points | update_JxW_values;
	const dealii::UpdateFlags faceUpdateFlags = update_values
			| update_quadrature_points | update_JxW_values
			| update_normal_vectors;
	const dealii::UpdateFlags neighborFaceUpdateFlags = update_values;
// Finite Element
	dealii::FEValues<dim> feCellValues(m_mapping, *m_fe, *m_quadrature,
			cellUpdateFlags);
	dealii::FEFaceValues<dim> feFaceValues(m_mapping, *m_fe, *m_faceQuadrature,
			faceUpdateFlags);
	dealii::FESubfaceValues<dim> feSubfaceValues(m_mapping, *m_fe,
			*m_faceQuadrature, faceUpdateFlags);
	dealii::FEFaceValues<dim> feNeighborFaceValues(m_mapping, *m_fe,
			*m_faceQuadrature, neighborFaceUpdateFlags);

	const size_t dofs_per_cell = m_fe->dofs_per_cell;
	const size_t n_quadrature_points = m_quadrature->size();

// Initialize matrices
	vector<double> localMassMatrix(dofs_per_cell);
	vector<dealii::FullMatrix<double> > localDerivativeMatrices;
	for (size_t i = 0; i < dim; i++) {
		dealii::FullMatrix<double> D_i(dofs_per_cell, dofs_per_cell);
		localDerivativeMatrices.push_back(D_i);
	}
	dealii::FullMatrix<double> localFaceMatrix(dofs_per_cell, dofs_per_cell);
	dealii::FullMatrix<double> localSystemMatrix(dofs_per_cell, dofs_per_cell);
	std::vector<types::global_dof_index> localDoFIndices(dofs_per_cell);

///////////////
// MAIN LOOP //
///////////////
	typename DoFHandler<dim>::active_cell_iterator cell =
			m_doFHandler->begin_active(), endc = m_doFHandler->end();
	for (; cell != endc; ++cell) {
		// calculate the fe values for the cell
		feCellValues.reinit(cell);

		// get global degrees of freedom
		cell->get_dof_indices(localDoFIndices);

		// assemble local cell matrices
		assembleLocalMassMatrix(feCellValues, dofs_per_cell,
				n_quadrature_points, localMassMatrix, localDoFIndices);
		for (size_t i = 0; i < dim; i++) {
			assembleLocalDerivativeMatrix(i, feCellValues, dofs_per_cell,
					n_quadrature_points, localDerivativeMatrices.at(i));
		}

		// assemble faces and put together
		for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
			// calculate local diagonal block (cell) matrix -D
			calculateAndDistributeLocalStiffnessMatrix(i, localDerivativeMatrices,
					localSystemMatrix, localDoFIndices, dofs_per_cell);
			// calculate face contributions M^-1 * R
			assembleAndDistributeLocalFaceMatrices(i,
					cell, feFaceValues, feSubfaceValues, feNeighborFaceValues,
					dofs_per_cell, n_quadrature_points, localFaceMatrix);
		}
	}
	// Mulitply by inverse mass matrix
	for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
		divideByDiagonalMassMatrix(m_systemMatrix.at(i), m_massMatrix);
	}
}
/// The template parameter must be made explicit in order for the code to compile
template void DataMinLee2011<2>::reassemble();
template void DataMinLee2011<3>::reassemble();

} /* namespace natrium */
