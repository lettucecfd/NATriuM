/**
 * @file SEDGMinLee.cpp
 * @short Advection scheme proposed by Min and Lee (2011)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SEDGMinLee.h"

#include "fstream"

#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/dofs/dof_renumbering.h>
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"

#include "../problemdescription/PeriodicBoundary.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
SEDGMinLee<dim>::SEDGMinLee(shared_ptr<Triangulation<dim> > triangulation,
		shared_ptr<BoundaryCollection<dim> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel,
		string inputDirectory, bool useCentralFlux) :
		m_tria(triangulation), m_boundaries(boundaries), m_mapping(), m_boltzmannModel(
				boltzmannModel), m_useCentralFlux(useCentralFlux) {
	// assertions
	assert(orderOfFiniteElement >= 2);

	// make dof handler
	m_quadrature = make_shared<QGaussLobatto<dim> >(orderOfFiniteElement);
	m_faceQuadrature = make_shared<QGaussLobatto<dim - 1> >(
			orderOfFiniteElement);
	m_fe = make_shared<FE_DGQArbitraryNodes<dim> >(
			QGaussLobatto<1>(orderOfFiniteElement));
	m_doFHandler = make_shared<DoFHandler<dim> >(*triangulation);

	for (size_t i = 0; i < m_boltzmannModel->getQ(); i++) {
		// (the first one will be empty all the time; just for consistency of indices with the distribution functions)
		m_systemMatrix.push_back(distributed_sparse_matrix());
	}
	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);
	updateSparsityPattern();

	// define relation between dofs and quadrature nodes
	m_facedof_to_q_index = map_facedofs_to_q_index();
	m_celldof_to_q_index = map_celldofs_to_q_index();
	m_q_index_to_facedof = map_q_index_to_facedofs();

	// reassemble or read file
	if (inputDirectory.empty()) {
		reassemble();
	} else {
		loadCheckpoint(inputDirectory);
	}

} /* SEDGMinLee<dim>::SEDGMinLee */
/// The template parameter must be made explicit in order for the code to compile
template SEDGMinLee<2>::SEDGMinLee(shared_ptr<Triangulation<2> > triangulation,
		shared_ptr<BoundaryCollection<2> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel,
		string inputDirectory, bool useCentralFlux);
template SEDGMinLee<3>::SEDGMinLee(shared_ptr<Triangulation<3> > triangulation,
		shared_ptr<BoundaryCollection<3> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<BoltzmannModel> boltzmannModel,
		string inputDirectory, bool useCentralFlux);



template<size_t dim>
void SEDGMinLee<dim>::reassemble() {
// TODO: if Triangulation changed: reinit dof-handler and sparsity pattern in some way

/////////////////////////////////
// Initialize Finite Element ////
/////////////////////////////////
// Define update flags (which values have to be known at each cell, face, neighbor face)
	const dealii::UpdateFlags cellUpdateFlags = update_values | update_gradients

	| update_quadrature_points | update_JxW_values | update_inverse_jacobians;
	const dealii::UpdateFlags faceUpdateFlags = update_values

	| update_quadrature_points | update_JxW_values | update_normal_vectors;
	const dealii::UpdateFlags neighborFaceUpdateFlags = update_values
			| update_JxW_values | update_normal_vectors;
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
		assembleLocalDerivativeMatrices(feCellValues, dofs_per_cell,
				n_quadrature_points, localDerivativeMatrices);

		// assemble faces and put together
		for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
// calculate local diagonal block (cell) matrix -D
			calculateAndDistributeLocalStiffnessMatrix(i,
					localDerivativeMatrices, localSystemMatrix, localDoFIndices,
					dofs_per_cell);
// calculate face contributions  R
			assembleAndDistributeLocalFaceMatrices(i, cell, feFaceValues,
					feSubfaceValues, feNeighborFaceValues, dofs_per_cell,
					n_quadrature_points, localFaceMatrix);
		}
	}
	// Mulitply by inverse mass matrix
	for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
		divideByDiagonalMassMatrix(m_systemMatrix.at(i), m_massMatrix);
	}
} /* reassemble */
/// The template parameter must be made explicit in order for the code to compile
template void SEDGMinLee<2>::reassemble();
template void SEDGMinLee<3>::reassemble();


template<size_t dim>
void SEDGMinLee<dim>::updateSparsityPattern() {

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
template void SEDGMinLee<2>::updateSparsityPattern();
//TODO generalize to 3D
//template void SEDGMinLee<3>::updateSparsityPattern();

template<size_t dim>
void SEDGMinLee<dim>::assembleLocalMassMatrix(
		const dealii::FEValues<dim>& feValues, size_t dofs_per_cell,
		size_t n_q_points, vector<double> &massMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs) {
	// initialize with zeros
	// TODO the fill operation can be cut out in the final implementation,
	// changing the += in the loop to =
	std::fill(massMatrix.begin(), massMatrix.end(), 0.0);

	// fill diagonal "matrix"
	for (size_t i = 0; i < dofs_per_cell; i++) {
		size_t q_point = m_celldof_to_q_index.at(i);
		massMatrix.at(i) += feValues.shape_value(i, q_point)
				* feValues.shape_value(i, q_point) * feValues.JxW(q_point);
	}

	// Assemble to global mass matrix
	for (size_t i = 0; i < dofs_per_cell; i++) {
		m_massMatrix(globalDoFs.at(i)) = massMatrix.at(i);
	}
} /*assembleLocalMassMatrix*/
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::assembleLocalMassMatrix(
		const dealii::FEValues<2>& feValues, size_t dofs_per_cell,
		size_t n_q_points, vector<double> &massMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs);
template void SEDGMinLee<3>::assembleLocalMassMatrix(
		const dealii::FEValues<3>& feValues, size_t dofs_per_cell,
		size_t n_q_points, vector<double> &massMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs);

template<size_t dim>
void SEDGMinLee<dim>::assembleLocalDerivativeMatrices(
		const dealii::FEValues<dim>& feValues, size_t dofs_per_cell,
		size_t n_q_points,
		vector<dealii::FullMatrix<double> >&derivativeMatrix) const {
	for (size_t i = 0; i < dim; i++) {
		derivativeMatrix.at(i) = 0;
	}
	for (size_t i = 0; i < dofs_per_cell; i++) {
		for (size_t j = 0; j < dofs_per_cell; j++) {
			// the shape value is zero for all q, except i<->q
			size_t q_point = m_celldof_to_q_index.at(i);
			Tensor<1, dim> integrandAtQ;
			integrandAtQ = feValues.shape_grad(j, q_point);
			integrandAtQ *= (feValues.shape_value(i, q_point)
					* feValues.JxW(q_point));
			for (size_t k = 0; k < dim; k++) {
				derivativeMatrix.at(k)(i, j) += integrandAtQ[k];
			}

		}
	}

} /* assembleLocalDerivativeMatrix */
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::assembleLocalDerivativeMatrices(
		const dealii::FEValues<2>& feValues, size_t dofs_per_cell,
		size_t n_q_points,
		vector<dealii::FullMatrix<double> > &derivativeMatrix) const;
template void SEDGMinLee<3>::assembleLocalDerivativeMatrices(
		const dealii::FEValues<3>& feValues, size_t dofs_per_cell,
		size_t n_q_points,
		vector<dealii::FullMatrix<double> > &derivativeMatrix) const;

template<size_t dim>
void SEDGMinLee<dim>::assembleAndDistributeLocalFaceMatrices(size_t i,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		dealii::FEFaceValues<dim>& feFaceValues,
		dealii::FESubfaceValues<dim>& feSubfaceValues,
		dealii::FEFaceValues<dim>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix) {
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
					assembleAndDistributeInternalFace(i, cell, j, neighborCell,
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
			assembleAndDistributeInternalFace(i, cell, j, neighbor,
					dealii::GeometryInfo<dim>::opposite_face[j], feFaceValues,
					feSubfaceValues, feNeighborFaceValues);
		} /* if (face at boundary) {} else {} */

	}

} /* assembleLocalFaceMatrix */
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::assembleAndDistributeLocalFaceMatrices(size_t i,
		typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		dealii::FEFaceValues<2>& feFaceValues,
		dealii::FESubfaceValues<2>& feSubfaceValues,
		dealii::FEFaceValues<2>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix);
template void SEDGMinLee<3>::assembleAndDistributeLocalFaceMatrices(size_t i,
		typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		dealii::FEFaceValues<3>& feFaceValues,
		dealii::FESubfaceValues<3>& feSubfaceValues,
		dealii::FEFaceValues<3>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix);

template<size_t dim>
void SEDGMinLee<dim>::divideByDiagonalMassMatrix(
		distributed_sparse_matrix& matrix,
		const distributed_vector& massMatrix) {
	size_t n = massMatrix.size();
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
// TODO Parallel loop
			if (m_sparsityPattern.exists(i, j)) {
				matrix.set(i, j, matrix(i, j) / massMatrix(i));
			}
		}
	}
}
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::divideByDiagonalMassMatrix(
		distributed_sparse_matrix& matrix,
		const distributed_vector& massMatrix);
template void SEDGMinLee<3>::divideByDiagonalMassMatrix(
		distributed_sparse_matrix& matrix,
		const distributed_vector& massMatrix);

template<> void SEDGMinLee<2>::calculateAndDistributeLocalStiffnessMatrix(
		size_t i, const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs,
		size_t dofsPerCell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	systemMatrix = derivativeMatrices.at(0);
	systemMatrix *= (-(m_boltzmannModel->getDirection(i)[0]));
	systemMatrix.add(-(m_boltzmannModel->getDirection(i)[1]),
			derivativeMatrices.at(1));
// distribute to global system matrix
	for (unsigned int j = 0; j < dofsPerCell; j++)
		for (unsigned int k = 0; k < dofsPerCell; k++) {
			m_systemMatrix.at(i).add(globalDoFs[j], globalDoFs[k],
					systemMatrix(j, k));
		}
}
template<> void SEDGMinLee<3>::calculateAndDistributeLocalStiffnessMatrix(
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
void SEDGMinLee<dim>::assembleAndDistributeInternalFace(size_t direction,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<dim>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		dealii::FESubfaceValues<dim>& feSubfaceValues,
		dealii::FEFaceValues<dim>& feNeighborFaceValues) {
	// get the required FE Values for the local cell
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<Point<dim> > &normals = feFaceValues.get_normal_vectors();

	// get the required dofs of the neighbor cell
	//	typename DoFHandler<dim>::face_iterator neighborFace = neighborCell->face(
	//			neighborFaceNumber);
	feNeighborFaceValues.reinit(neighborCell, neighborFaceNumber);

	// TODO clean up construction; or (better): assembly directly
	dealii::FullMatrix<double> cellFaceMatrix(feFaceValues.dofs_per_cell);
	dealii::FullMatrix<double> neighborFaceMatrix(
			feNeighborFaceValues.dofs_per_cell, feFaceValues.dofs_per_cell);
	cellFaceMatrix = 0;
	neighborFaceMatrix = 0;

	// loop over all quadrature points at the face
	for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
		size_t thisDoF = m_q_index_to_facedof.at(faceNumber).at(q);
		size_t neighborDoF = m_q_index_to_facedof.at(neighborFaceNumber).at(q);

		// calculate matrix entries
		vector<double> factor(JxW);
		double exn = 0.0;
		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
			exn += normals.at(q)(i)
					* m_boltzmannModel->getDirection(direction)(i);
		}
		for (size_t i = 0; i < factor.size(); i++) {
			factor.at(i) *= exn;
		}

		if (m_useCentralFlux) {

			cellFaceMatrix(thisDoF, thisDoF) = 0.5 * factor.at(q);
			neighborFaceMatrix(thisDoF, neighborDoF) = 0.5 * factor.at(q);

		} else if (exn < 0) { // otherwise: no contributions

			cellFaceMatrix(thisDoF, thisDoF) = factor.at(q);
			neighborFaceMatrix(thisDoF, neighborDoF) = -factor.at(q);
		}
	}

	// get DoF indices
	// TODO cut out construction (allocation); Allocating two vectors in most inner loop is too expensive
	vector<dealii::types::global_dof_index> localDoFIndices(
			feFaceValues.dofs_per_cell);
	vector<dealii::types::global_dof_index> neighborDoFIndices(
			feFaceValues.dofs_per_cell);
	cell->get_dof_indices(localDoFIndices);
	neighborCell->get_dof_indices(neighborDoFIndices);

	/// Distribute to global matrix
	for (size_t i = 0; i < feFaceValues.dofs_per_cell; i++) {
		for (size_t j = 0; j < feFaceValues.dofs_per_cell; j++) {
			m_systemMatrix.at(direction).add(localDoFIndices[i],
					localDoFIndices[j], cellFaceMatrix(i, j));
			m_systemMatrix.at(direction).add(localDoFIndices[i],
					neighborDoFIndices[j], neighborFaceMatrix(i, j));
		}
	}

// WARNING: Hanging nodes are not implemented, yet
// TODO Implement local refinement

} /* assembleAndDistributeInternalFace */
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::assembleAndDistributeInternalFace(size_t direction,
		typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<2>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValues<2>& feFaceValues,
		dealii::FESubfaceValues<2>& feSubfaceValues,
		dealii::FEFaceValues<2>& feNeighborFaceValues);
template void SEDGMinLee<3>::assembleAndDistributeInternalFace(size_t direction,
		typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<3>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValues<3>& feFaceValues,
		dealii::FESubfaceValues<3>& feSubfaceValues,
		dealii::FEFaceValues<3>& feNeighborFaceValues);

template<size_t dim>
std::map<size_t, size_t> natrium::SEDGMinLee<dim>::map_celldofs_to_q_index() const {
	const dealii::UpdateFlags cellUpdateFlags = update_values
			| update_quadrature_points;
	// Finite Element
	dealii::FEValues<dim> feCellValues(m_mapping, *m_fe, *m_quadrature,
			cellUpdateFlags);
	const size_t dofs_per_cell = m_fe->dofs_per_cell;
	const size_t n_quadrature_points = m_quadrature->size();
	// take first cell
	typename DoFHandler<dim>::active_cell_iterator cell =
			m_doFHandler->begin_active();
	std::map<size_t, size_t> result;

	/// find quadrature node for every DoF
	feCellValues.reinit(cell);
	for (size_t i = 0; i < dofs_per_cell; i++) {
		int unique = 0;
		for (size_t q = 0; q < n_quadrature_points; q++) {
			if (feCellValues.shape_value(i, q) > 1e-10) {
				assert(fabs(feCellValues.shape_value(i, q) - 1) < 1e-10);
				unique += 1;
				result.insert(std::make_pair(i, q));
			}
		}
		assert(unique == 1);
	}
	return result;
}
// The template parameter must be made explicit in order for the code to compile.
template std::map<size_t, size_t> SEDGMinLee<2>::map_celldofs_to_q_index() const;
template std::map<size_t, size_t> SEDGMinLee<3>::map_celldofs_to_q_index() const;

template<size_t dim>
vector<std::map<size_t, size_t> > natrium::SEDGMinLee<dim>::map_facedofs_to_q_index() const {
	const dealii::UpdateFlags faceUpdateFlags = update_values
			| update_quadrature_points;
	dealii::FEFaceValues<dim> feFaceValues(m_mapping, *m_fe, *m_faceQuadrature,
			faceUpdateFlags);

	typename DoFHandler<dim>::active_cell_iterator cell =
			m_doFHandler->begin_active();
	// LOOP over all faces
	vector<std::map<size_t, size_t> > result;
	for (size_t f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
		feFaceValues.reinit(cell, f);
		std::map<size_t, size_t> resultForFaceF;
		for (size_t i = 0; i < m_fe->dofs_per_cell; i++) {
			int unique = 0;
			for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
				if (feFaceValues.shape_value(i, q) > 1e-10) {
					assert(fabs(feFaceValues.shape_value(i, q) - 1) < 1e-10);
					unique += 1;
					resultForFaceF.insert(std::make_pair(i, q));
				}
			}
// Test, if the relationship doF <-> quadrature points in unique
			assert(unique <= 1);
		}
		result.push_back(resultForFaceF);
		assert(resultForFaceF.size() == feFaceValues.n_quadrature_points);
	}
	return result;
}
// The template parameter must be made explicit in order for the code to compile.
template vector<std::map<size_t, size_t> > SEDGMinLee<2>::map_facedofs_to_q_index() const;
template vector<std::map<size_t, size_t> > SEDGMinLee<3>::map_facedofs_to_q_index() const;

template<size_t dim>
vector<std::map<size_t, size_t> > natrium::SEDGMinLee<dim>::map_q_index_to_facedofs() const {
	const dealii::UpdateFlags faceUpdateFlags = update_values
			| update_quadrature_points;
	dealii::FEFaceValues<dim> feFaceValues(m_mapping, *m_fe, *m_faceQuadrature,
			faceUpdateFlags);

	typename DoFHandler<dim>::active_cell_iterator cell =
			m_doFHandler->begin_active();
	// LOOP over all faces
	vector<std::map<size_t, size_t> > result;
	for (size_t f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
		feFaceValues.reinit(cell, f);
		std::map<size_t, size_t> resultForFaceF;
		for (size_t i = 0; i < m_fe->dofs_per_cell; i++) {
			int unique = 0;
			for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
				if (feFaceValues.shape_value(i, q) > 1e-10) {
					assert(fabs(feFaceValues.shape_value(i, q) - 1) < 1e-10);
					unique += 1;
					resultForFaceF.insert(std::make_pair(q, i));
				}
			}
// Test, if the relationship doF <-> quadrature points in unique
			assert(unique <= 1);
		}
		result.push_back(resultForFaceF);
		assert(resultForFaceF.size() == feFaceValues.n_quadrature_points);
	}
	return result;
} /* map_q_index_to_facedofs */
// The template parameter must be made explicit in order for the code to compile.
template vector<std::map<size_t, size_t> > SEDGMinLee<2>::map_q_index_to_facedofs() const;
template vector<std::map<size_t, size_t> > SEDGMinLee<3>::map_q_index_to_facedofs() const;


template<size_t dim>
void SEDGMinLee<dim>::stream() {
}
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::stream();
template void SEDGMinLee<3>::stream();


template<size_t dim>
void SEDGMinLee<dim>::saveMatricesToFiles(const string& directory) const {
	// write system matrices to files
	try {
		for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
			// filename
			std::stringstream filename;
			filename << directory << "/checkpoint_system_matrix_" << i << ".dat";
			std::ofstream file(filename.str().c_str());
			m_systemMatrix.at(i).block_write(file);
		}
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while writing the system matrices to files: Please make shure you have writing permission. Quick fix: Remove StreamingMatrices from OutputFlags");
	}
	// Write the mass matrix
	try {
		// filename
		std::stringstream filename;
		filename << directory << "/checkpoint_mass_matrix.dat";
		std::ofstream file(filename.str().c_str());
		m_massMatrix.block_write(file);
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while writing the mass matrix to file: Please make shure you have writing permission. Quick fix: Remove StreamingMatrices from OutputFlags");
	}
}
template void SEDGMinLee<2>::saveMatricesToFiles(const string& directory) const;
template void SEDGMinLee<3>::saveMatricesToFiles(const string& directory) const;

template<size_t dim>
void SEDGMinLee<dim>::loadMatricesFromFiles(const string& directory) {
	// read the system matrices from file
	try {
		for (size_t i = 1; i < m_boltzmannModel->getQ(); i++) {
			// filename
			std::stringstream filename;
			filename << directory << "/checkpoint_system_matrix_" << i << ".dat";
			std::ifstream file(filename.str().c_str());
			m_systemMatrix.at(i).block_read(file);
		}
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while reading the system matrices from file: Please switch off the restart option to start the simulation from the beginning.");
	}
	// Read the mass matrix
	try {
		// filename
		std::stringstream filename;
		filename << directory << "/checkpoint_mass_matrix.dat";
		std::ifstream file(filename.str().c_str());
		m_massMatrix.block_read(file);
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while reading the mass matrix from file: Please switch off the restart option to start the simulation from the beginning.");
	}
	// TODO Test: Is the matrix OK?
}
template void SEDGMinLee<2>::loadMatricesFromFiles(const string& directory);
template void SEDGMinLee<3>::loadMatricesFromFiles(const string& directory);

template<size_t dim> void SEDGMinLee<dim>::writeStatus(
		const string& directory) const {
	//make file
	std::stringstream filename;
	filename << directory << "/checkpoint_status.dat";
	std::ofstream outfile(filename.str().c_str());

	//write number of cells
	outfile << m_tria->n_cells() << endl;
	//write order of fe
	outfile << m_fe->get_degree() << endl;
	//write number of dofs
	outfile << this->getNumberOfDoFs() << endl;
	//write D
	outfile << m_boltzmannModel->getD() << endl;
	//write Q
	outfile << m_boltzmannModel->getQ() << endl;
	//write magic number of cell geometry
	outfile << calcMagicNumber() << endl;
	//write dqScaling1
	outfile << m_boltzmannModel->getDirection(1)(0) << endl;
	//write dqScaling2
	outfile << m_boltzmannModel->getDirection(1)(1) << endl;
	//write fluxType
	outfile << m_useCentralFlux << endl;
	//write advectionType
	outfile << "SEDGMinLee" << endl;
}
template void SEDGMinLee<2>::writeStatus(const string& directory) const;
template void SEDGMinLee<3>::writeStatus(const string& directory) const;

template<size_t dim> bool SEDGMinLee<dim>::isStatusOK(const string& directory,
		string& message) const {
	//read file
	std::stringstream filename;
	filename << directory << "/checkpoint_status.dat";
	std::ifstream infile(filename.str().c_str());

	// check if status file exists
	if (not infile) {
		message = "No checkpoint found. Please disable restart option.";
		return false;
	}
	//number of cells
	size_t tmp;
	infile >> tmp;
	if (tmp != m_tria->n_cells()) {
		message = "Number of cells not equal.";
		return false;
	}
	// order of fe
	infile >> tmp;
	if (tmp != m_fe->get_degree()) {
		message = "Order of finite element not equal.";
		return false;
	}
	// number of dofs
	infile >> tmp;
	if (tmp != this->getNumberOfDoFs()) {
		message = "Number of degrees of freedom not equal.";
		return false;
	}
	// D
	infile >> tmp;
	if (tmp != m_boltzmannModel->getD()) {
		message = "Dimension not equal.";
		return false;
	}
	// Q
	infile >> tmp;
	if (tmp != m_boltzmannModel->getQ()) {
		message = "Number of particle velocities not equal.";
		return false;
	}
	// magic number of cell geometry
	double dtmp;
	infile >> dtmp;
	if (fabs(dtmp - calcMagicNumber()) > 1e-1) {
		message = "Triangulation (or at least its magic number) not equal.";
		return false;
	}
	// dqScaling1
	infile >> dtmp;
	if (fabs(dtmp - m_boltzmannModel->getDirection(1)(0)) > 1e-5) {
		message = "Scaling of Boltzmann model (1st coordinate) not equal.";
		return false;
	}
	// dqScaling2
	infile >> dtmp;
	if (fabs(dtmp - m_boltzmannModel->getDirection(1)(1)) > 1e-5) {
		message = "Scaling of Boltzmann model (2nd) not equal.";
		return false;
	}
	// fluxType
	infile >> tmp;
	if (tmp != m_useCentralFlux) {
		message = "Flux not equal.";
		return false;
	}
	// advectionType
	string stmp;
	infile >> stmp;
	if (stmp != "SEDGMinLee") {
		message = "AdvectionOperator Type not equal.";
		return false;
	}

	return true;
}
template bool SEDGMinLee<2>::isStatusOK(const string& directory,
		string& message) const;
template bool SEDGMinLee<3>::isStatusOK(const string& directory,
		string& message) const;

template<size_t dim>
double SEDGMinLee<dim>::calcMagicNumber() const {

	double magicNumber = m_fe->get_degree() / m_tria->n_cells();

	vector<dealii::Point<dim> > supportPoints(this->getNumberOfDoFs());
	mapDoFsToSupportPoints(supportPoints);
	for (size_t i = 0; i < supportPoints.size();
			i += std::max(1., supportPoints.size() / 100.)) {
		magicNumber += (1 + (i % 1000)) / (1 + i) * 1.
				/ (supportPoints.at(i)(0) + 1) * 1.
				/ (1 + supportPoints.at(i)(1)) * 1.
				/ (1 + supportPoints.at(i)(dim - 1));
	}
	return magicNumber;
}

template double SEDGMinLee<2>::calcMagicNumber() const;
template double SEDGMinLee<3>::calcMagicNumber() const;

template<size_t dim>
void SEDGMinLee<dim>::loadCheckpoint(const string& directory){
	// check if stuff can be read from file. Else throw exception
	string message;
	if (isStatusOK(directory, message)) {
		loadMatricesFromFiles(directory);
	} else {
		std::stringstream errorMessage;
		errorMessage << "Restart not possible. " << message;
		throw AdvectionSolverException(errorMessage.str().c_str());
	}
}
template void SEDGMinLee<2>::loadCheckpoint(const string& directory);
template void SEDGMinLee<3>::loadCheckpoint(const string& directory);


template<size_t dim>
void SEDGMinLee<dim>::saveCheckpoint(const string& directory) const{
	writeStatus(directory);
	saveMatricesToFiles(directory);
}
template void SEDGMinLee<2>::saveCheckpoint(const string& directory) const;
template void SEDGMinLee<3>::saveCheckpoint(const string& directory) const;

} /* namespace natrium */
