/**
 * @file SEDGMinLee.cpp
 * @short Advection scheme proposed by Min and Lee (2011)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SEDGMinLee.h"

#include "fstream"

#include "deal.II/lac/compressed_sparsity_pattern.h"
#include "deal.II/dofs/dof_renumbering.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/lac/matrix_iterator.h"

#include "../problemdescription/PeriodicBoundary.h"
#include "../problemdescription/MinLeeBoundary.h"

#include "../stencils/Stencil.h"

#include "../utilities/DealiiExtensions.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
SEDGMinLee<dim>::SEDGMinLee(shared_ptr<Mesh<dim> > triangulation,
		shared_ptr<BoundaryCollection<dim> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<Stencil> Stencil,
		string inputDirectory, bool useCentralFlux) :
		m_tria(triangulation), m_boundaries(boundaries), m_mapping(orderOfFiniteElement), m_stencil(
				Stencil), m_useCentralFlux(useCentralFlux) {
	// assertions
	assert(orderOfFiniteElement >= 1);
	assert(Stencil->getD() == dim);

	// make dof handler
	m_quadrature = make_shared<QGaussLobatto<dim> >(orderOfFiniteElement + 1);
	m_faceQuadrature = make_shared<QGaussLobatto<dim - 1> >(
			orderOfFiniteElement + 1);
	m_fe = make_shared<FE_DGQArbitraryNodes<dim> >(
			QGaussLobatto<1>(orderOfFiniteElement + 1));
	m_doFHandler = make_shared<DoFHandler<dim> >(*triangulation);

	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);

#ifdef WITH_TRILINOS_MPI
	//get locally owned and locally relevant dofs
	m_locallyOwnedDofs = m_doFHandler->locally_owned_dofs ();
	DoFTools::extract_locally_relevant_dofs (*m_doFHandler,
			m_locallyRelevantDofs);
#endif
	updateSparsityPattern();

	// define relation between dofs and quadrature nodes
	m_facedof_to_q_index = map_facedofs_to_q_index();
	m_celldof_to_q_index = map_celldofs_to_q_index();
	m_q_index_to_facedof = map_q_index_to_facedofs();

	// set size for the system vector
#ifdef WITH_TRILINOS
	m_systemVector.reinit(m_stencil->getQ() - 1);
	for (size_t i = 0; i < m_stencil->getQ() - 1; i++) {
#ifdef WITH_TRILINOS_MPI
		m_systemVector.block(i).reinit(m_locallyOwnedDofs, m_locallyRelevantDofs, MPI_COMM_WORLD);
#else
		m_systemVector.block(i).reinit(m_doFHandler->n_dofs());
#endif
	}
	m_systemVector.collect_sizes();
#else
	m_systemVector.reinit(m_stencil->getQ() - 1, m_doFHandler->n_dofs());
#endif
	// reassemble or read file
	if (inputDirectory.empty()) {
		reassemble();
	} else {
		loadCheckpoint(inputDirectory);
	}

} /* SEDGMinLee<dim>::SEDGMinLee */
/// The template parameter must be made explicit in order for the code to compile
template SEDGMinLee<2>::SEDGMinLee(shared_ptr<Mesh<2> > triangulation,
		shared_ptr<BoundaryCollection<2> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<Stencil> Stencil,
		string inputDirectory, bool useCentralFlux);
template SEDGMinLee<3>::SEDGMinLee(shared_ptr<Mesh<3> > triangulation,
		shared_ptr<BoundaryCollection<3> > boundaries,
		size_t orderOfFiniteElement, shared_ptr<Stencil> Stencil,
		string inputDirectory, bool useCentralFlux);

template<size_t dim>
void SEDGMinLee<dim>::reassemble() {
// TODO: if Mesh changed: reinit dof-handler and sparsity pattern in some way

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
	vector<double> inverseLocalMassMatrix(dofs_per_cell);
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

		// invert local mass matrix
		for (size_t i = 0; i < dofs_per_cell; i++){
			inverseLocalMassMatrix.at(i) = 1./localMassMatrix.at(i);
		}

		// assemble faces and put together
		for (size_t alpha = 1; alpha < m_stencil->getQ(); alpha++) {
// calculate local diagonal block (cell) matrix -D
			calculateAndDistributeLocalStiffnessMatrix(alpha,
					localDerivativeMatrices, localSystemMatrix, inverseLocalMassMatrix, localDoFIndices,
					dofs_per_cell);
// calculate face contributions  R
			assembleAndDistributeLocalFaceMatrices(alpha, cell, feFaceValues,
					feSubfaceValues, feNeighborFaceValues, dofs_per_cell,
					n_quadrature_points, localFaceMatrix, inverseLocalMassMatrix);
		}
	}

//#endif

} /* reassemble */
/// The template parameter must be made explicit in order for the code to compile
template void SEDGMinLee<2>::reassemble();
template void SEDGMinLee<3>::reassemble();

template<size_t dim>
void SEDGMinLee<dim>::updateSparsityPattern() {
	// TODO only update sparsity pattern for changed cells
	///////////////////////////////////////////////////////
	// Setup sparsity pattern (completely manually):
	///////////////////////////////////////////////////////
	// allocate sizes
	size_t n_blocks = m_stencil->getQ() - 1;
	size_t n_dofs_per_block = m_doFHandler->n_dofs();
	
	const dealii::UpdateFlags faceUpdateFlags = update_values
			| update_quadrature_points;
	dealii::FEFaceValues<dim>* feFaceValues = new FEFaceValues<dim>(m_mapping,
			*m_fe, *m_faceQuadrature, faceUpdateFlags);

#ifdef WITH_TRILINOS
	// Trilinos can work with improved sparsity structures

	CompressedSparsityPattern cSparseDiag(n_dofs_per_block, n_dofs_per_block);
	CompressedSparsityPattern cSparseOpposite(n_dofs_per_block, n_dofs_per_block);
	CompressedSparsityPattern cSparseEmpty(n_dofs_per_block, n_dofs_per_block);



	// create cell maps for periodic boundary
	for (typename BoundaryCollection<dim>::ConstPeriodicIterator periodic =
			m_boundaries->getPeriodicBoundaries().begin();
			periodic != m_boundaries->getPeriodicBoundaries().end();
			periodic++) {
		// Periodic boundaries have two boundary indicators; (and are stored twice in the map)
		// skip double execution of addToSparsityPattern
		if (periodic->first == periodic->second->getBoundaryIndicator1()) {
			periodic->second->createCellMap(*m_doFHandler);
		}
	}

	//reorder degrees of freedom
	//DoFRenumbering::Cuthill_McKee(*m_doFHandler);
	//The Renumbering operation is commented out because of its quadratic complexity.
	//When it is used, it takes the most time of all assembly functions for meshes with >250 cells
	//The main point of renumbering algorithms is to make solving linear equation systems faster.
	//As we have no LES to solve, renumbering does not bring anything here.
	//However, when porting the code to distributed triangulations, it might become an issue again,
	//as the numbering of the DoFs determines the partition (and thus the number of halo nodes
	//for each MPI process). Furthermore, using implicit time integration schemes could bring
	//the issue of renumbering up again, as they require linear equation systems to be solved.

	// make diagonal block 0,0 which can be copied to the other ones
	DealIIExtensions::make_sparser_flux_sparsity_pattern(*m_doFHandler,
			cSparseDiag, *m_boundaries, feFaceValues);
	delete feFaceValues;
	/*DoFTools::make_flux_sparsity_pattern(*m_doFHandler,
	 cSparse.block(0, 0));*/

	// add entries for non-periodic boundaries
	for (typename BoundaryCollection<dim>::ConstMinLeeIterator minLeeIterator =
			m_boundaries->getMinLeeBoundaries().begin();
			minLeeIterator != m_boundaries->getMinLeeBoundaries().end();
			minLeeIterator++) {
		minLeeIterator->second->addToSparsityPattern(cSparseOpposite,
				*m_doFHandler, *m_stencil);
	}

	//reinitialize matrices
	//In order to store the sparsity pattern for blocks with same pattern only once: initialize from other block
	m_systemMatrix.reinit(n_blocks, n_blocks);
	m_systemMatrix.block(0, 0).reinit(cSparseDiag);
	size_t first_opposite = m_stencil->getIndexOfOppositeDirection(1) - 1;
	size_t some_empty = m_stencil->getIndexOfOppositeDirection(1);
	m_systemMatrix.block(0, some_empty).reinit(cSparseEmpty);
	m_systemMatrix.block(0, first_opposite).reinit(cSparseOpposite);
	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			if ((I == 0) and (J == 0)) {
				continue;
			}
			if ((I == 0) and (J == some_empty)) {
				continue;
			}
			if ((I == 0) and (J == first_opposite)) {
				continue;
			}
			if (I == J) {
				m_systemMatrix.block(I, J).reinit(m_systemMatrix.block(0, 0));
				continue;
			}
			if (I == m_stencil->getIndexOfOppositeDirection(J + 1) - 1) {
				m_systemMatrix.block(I, J).reinit(m_systemMatrix.block(0, first_opposite));
				continue;
			} else {
				m_systemMatrix.block(I,J).reinit(m_systemMatrix.block(0,some_empty));
			}

		}
	}
	m_systemMatrix.collect_sizes();

#else
	BlockCompressedSparsityPattern cSparse(n_blocks, n_blocks);
	// TODO do not initialize empty blocks?
	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			cSparse.block(I, J).reinit(n_dofs_per_block, n_dofs_per_block);
		}
	}
	cSparse.collect_sizes();

	// THIS IS JUST A WORKAROUND
	// TODO get sparser sparsity pattern working also for simulations without trilinos.
	//DealIIExtensions::make_sparser_flux_sparsity_pattern(*m_doFHandler,
	//		cSparse.block(0, 0), *m_boundaries, feFaceValues);
	delete feFaceValues ;
	DoFTools::make_flux_sparsity_pattern(*m_doFHandler,
	 cSparse.block(0, 0));

	// add periodic boundaries to intermediate flux sparsity pattern
	size_t dofs_per_cell = m_doFHandler->get_fe().dofs_per_cell;
	for (typename BoundaryCollection<dim>::ConstPeriodicIterator periodic =
			m_boundaries->getPeriodicBoundaries().begin();
			periodic != m_boundaries->getPeriodicBoundaries().end();
			periodic++) {
		// Periodic boundaries have two boundary indicators; (and are stored twice in the map)
		// skip double execution of addToSparsityPattern
		if (periodic->first == periodic->second->getBoundaryIndicator1()) {
			periodic->second->createCellMap(*m_doFHandler);
			size_t only_on_block_0_0 = 1;
			periodic->second->addToSparsityPattern(cSparse, only_on_block_0_0,
					n_dofs_per_block, dofs_per_cell);
		}
	}

	// add entries for non-periodic boundaries
	/*for (typename BoundaryCollection<dim>::ConstMinLeeIterator minLeeIterator =
			m_boundaries->getMinLeeBoundaries().begin();
			minLeeIterator != m_boundaries->getMinLeeBoundaries().end();
			minLeeIterator++) {
		minLeeIterator->second->addToSparsityPattern(cSparse, *m_doFHandler,
				*m_stencil);
	}*/

	// initialize (static) sparsity pattern
	m_sparsityPattern.reinit(n_blocks, n_blocks);
	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			// the following distinction is valid because non-periodic boundaries
			// do not affect the diagonal blocks
			if ((I != J) or (I == 0)) {
				m_sparsityPattern.block(I, J).copy_from(cSparse.block(I, J));
			} else {
				m_sparsityPattern.block(I, J).copy_from(cSparse.block(0, 0));
			}
		}
	}
	m_sparsityPattern.collect_sizes();
	//reinitialize matrices
	m_systemMatrix.reinit(m_sparsityPattern);
#endif


}
/* updateSparsityPattern */
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
	std::fill(massMatrix.begin(), massMatrix.end(), 0.0);

// fill diagonal "matrix"
	for (size_t i = 0; i < dofs_per_cell; i++) {
		size_t q_point = m_celldof_to_q_index.at(i);
		massMatrix.at(i) += feValues.shape_value(i, q_point)
				* feValues.shape_value(i, q_point) * feValues.JxW(q_point);
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
void SEDGMinLee<dim>::assembleAndDistributeLocalFaceMatrices(size_t alpha,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		dealii::FEFaceValues<dim>& feFaceValues,
		dealii::FESubfaceValues<dim>& feSubfaceValues,
		dealii::FEFaceValues<dim>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix, const vector<double>& inverseLocalMassMatrix) {

// loop over all faces
	for (size_t j = 0; j < dealii::GeometryInfo<dim>::faces_per_cell; j++) {
		//Faces at boundary
		if (cell->face(j)->at_boundary()) {
			size_t boundaryIndicator = cell->face(j)->boundary_indicator();
			if (m_boundaries->isPeriodic(boundaryIndicator)) {
				// Apply periodic boundaries
				const shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
						m_boundaries->getPeriodicBoundary(boundaryIndicator);
				assert(periodicBoundary->isFaceInBoundary(cell, j));
				typename dealii::DoFHandler<dim>::cell_iterator neighborCell;
				size_t opposite_face =
						periodicBoundary->getOppositeCellAtPeriodicBoundary(
								cell, neighborCell);
				assembleAndDistributeInternalFace(alpha, cell, j, neighborCell,
						opposite_face, feFaceValues, feSubfaceValues,
						feNeighborFaceValues, inverseLocalMassMatrix);
			} else /* if is not periodic */{
				// Apply other boundaries
				if (typeid(*(m_boundaries->getBoundary(boundaryIndicator)))
						== typeid(MinLeeBoundary<dim> )) {
					const shared_ptr<MinLeeBoundary<dim> >& minLeeBoundary =
							m_boundaries->getMinLeeBoundary(boundaryIndicator);
					minLeeBoundary->assembleBoundary(alpha, cell, j,
							feFaceValues, *m_stencil,
							m_q_index_to_facedof.at(j), inverseLocalMassMatrix, m_systemMatrix,
							m_systemVector, m_useCentralFlux);
				}
			} /* endif isPeriodic */

		} else /* if is not at boundary */{
// Internal faces
			typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(
					j);
			assembleAndDistributeInternalFace(alpha, cell, j, neighbor,
					cell->neighbor_face_no(j), feFaceValues,
					feSubfaceValues, feNeighborFaceValues, inverseLocalMassMatrix);
		} /* endif is face at boundary*/

	}

} /* assembleLocalFaceMatrix */
// The template parameter must be made explicit in order for the code to compile.
template void SEDGMinLee<2>::assembleAndDistributeLocalFaceMatrices(
		size_t alpha,
		typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		dealii::FEFaceValues<2>& feFaceValues,
		dealii::FESubfaceValues<2>& feSubfaceValues,
		dealii::FEFaceValues<2>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix, const vector<double>& inverseLocalMassMatrix);
template void SEDGMinLee<3>::assembleAndDistributeLocalFaceMatrices(
		size_t alpha,
		typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		dealii::FEFaceValues<3>& feFaceValues,
		dealii::FESubfaceValues<3>& feSubfaceValues,
		dealii::FEFaceValues<3>& feNeighborFaceValues, size_t dofs_per_cell,
		size_t n_q_points, dealii::FullMatrix<double>& faceMatrix, const vector<double>& inverseLocalMassMatrix);

template<> void SEDGMinLee<2>::calculateAndDistributeLocalStiffnessMatrix(
		size_t alpha,
		const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix, const vector<double>& inverseLocalMassMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs,
		size_t dofsPerCell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	systemMatrix = derivativeMatrices.at(0);
	systemMatrix *= (-(m_stencil->getDirection(alpha)[0]));
	systemMatrix.add(-(m_stencil->getDirection(alpha)[1]),
			derivativeMatrices.at(1));
// distribute to global system matrix
	// avoid to call block() too often
	distributed_sparse_matrix& block = m_systemMatrix.block(alpha - 1,
			alpha - 1);
	for (unsigned int j = 0; j < dofsPerCell; j++)
		for (unsigned int k = 0; k < dofsPerCell; k++) {
			block.add(globalDoFs[j], globalDoFs[k], systemMatrix(j, k) * inverseLocalMassMatrix.at(j));
		}
}
template<> void SEDGMinLee<3>::calculateAndDistributeLocalStiffnessMatrix(
		size_t alpha,
		const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix, const vector<double>& inverseLocalMassMatrix,
		const std::vector<dealii::types::global_dof_index>& globalDoFs,
		size_t dofsPerCell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	systemMatrix = derivativeMatrices.at(0);
	systemMatrix *= (-m_stencil->getDirection(alpha)[0]);
	systemMatrix.add(-m_stencil->getDirection(alpha)[1],
			derivativeMatrices.at(1), -m_stencil->getDirection(alpha)[2],
			derivativeMatrices.at(2));
// distribute to global system matrix
	distributed_sparse_matrix& block = m_systemMatrix.block(alpha - 1,
			alpha - 1);
	for (unsigned int j = 0; j < dofsPerCell; j++)
		for (unsigned int k = 0; k < dofsPerCell; k++)
			block.add(globalDoFs[j], globalDoFs[k], systemMatrix(j, k) * inverseLocalMassMatrix.at(j));
}

template<size_t dim>
void SEDGMinLee<dim>::assembleAndDistributeInternalFace(size_t alpha,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<dim>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		dealii::FESubfaceValues<dim>& feSubfaceValues,
		dealii::FEFaceValues<dim>& feNeighborFaceValues, const vector<double>& inverseLocalMassMatrix) {
// get the required FE Values for the local cell
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<Point<dim> > &normals = feFaceValues.get_normal_vectors();

	if (4 == alpha){

	}
// get the required dofs of the neighbor cell
//	typename DoFHandler<dim>::face_iterator neighborFace = neighborCell->face(
//			neighborFaceNumber);
	feNeighborFaceValues.reinit(neighborCell, neighborFaceNumber);

	vector<dealii::types::global_dof_index> localDoFIndices(
			feFaceValues.get_fe().dofs_per_cell);
	vector<dealii::types::global_dof_index> neighborDoFIndices(
			feFaceValues.get_fe().dofs_per_cell);
	cell->get_dof_indices(localDoFIndices);
	neighborCell->get_dof_indices(neighborDoFIndices);

	// avoid many calls to block()
	distributed_sparse_matrix& block = m_systemMatrix.block(alpha - 1,
			alpha - 1);

// loop over all quadrature points at the face
	for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
		size_t thisDoF = m_q_index_to_facedof.at(faceNumber).at(q);
		assert (feFaceValues.shape_value(thisDoF, q) > 0);
		size_t neighborDoF = m_q_index_to_facedof.at(neighborFaceNumber).at(q);
		assert (feNeighborFaceValues.shape_value(neighborDoF, q) > 0);

		double cell_entry = 0.0;
		double neighbor_entry = 0.0;

		// calculate matrix entries
		double prefactor = JxW.at(q);
		double exn = 0.0;

		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {		// TODO efficient multiplication
			exn += normals.at(q)(i) * m_stencil->getDirection(alpha)(i);
		}
		prefactor *= exn;

		if (m_useCentralFlux) {
			cell_entry = 0.5 * prefactor;
			neighbor_entry = 0.5 * prefactor;

		} else if (exn < 0) { // otherwise: no contributions
			cell_entry = prefactor;
			neighbor_entry = -prefactor;
		}

		block.add(localDoFIndices[thisDoF], localDoFIndices[thisDoF],
				cell_entry * inverseLocalMassMatrix.at(thisDoF));
		block.add(localDoFIndices[thisDoF], neighborDoFIndices[neighborDoF],
				neighbor_entry * inverseLocalMassMatrix.at(thisDoF));
	}

// get DoF indices
// TODO cut out construction (allocation); Allocating two vectors in most inner loop is too expensive

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
		dealii::FEFaceValues<2>& feNeighborFaceValues, const vector<double>& inverseLocalMassMatrix);
template void SEDGMinLee<3>::assembleAndDistributeInternalFace(size_t direction,
		typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<3>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValues<3>& feFaceValues,
		dealii::FESubfaceValues<3>& feSubfaceValues,
		dealii::FEFaceValues<3>& feNeighborFaceValues, const vector<double>& inverseLocalMassMatrix);

template<size_t dim>
std::map<size_t, size_t> SEDGMinLee<dim>::map_celldofs_to_q_index() const {
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
vector<std::map<size_t, size_t> > SEDGMinLee<dim>::map_facedofs_to_q_index() const {
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
vector<std::map<size_t, size_t> > SEDGMinLee<dim>::map_q_index_to_facedofs() const {
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
		for (size_t i = 0; i < m_stencil->getQ() - 1; i++) {
			for (size_t j = 0; j < m_stencil->getQ() - 1; j++) {
				// filename
				std::stringstream filename;
				filename << directory << "/checkpoint_system_matrix_" << i
#ifdef WITH_TRILINOS_MPI
						<< "_" << j << "_" << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) << ".dat";
#else
						<< "_" << j << ".dat";
#endif
				std::ofstream file(filename.str().c_str());
#ifndef WITH_TRILINOS
				m_systemMatrix.block(i, j).block_write(file);
#else
				// TODO use trilinos functions for read and write. This here is really bad
#endif
			}
		}

	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while writing the system matrices to files: Please make sure you have writing permission. Quick fix: Remove StreamingMatrices from OutputFlags");
	}

// Write the system vector
	try {
		// filename
		std::stringstream filename;
		filename << directory << "/checkpoint_system_vector.dat";
		std::ofstream file(filename.str().c_str());
#ifdef WITH_TRILINOS
		for (size_t i = 0; i < m_stencil->getQ() - 1; i++) {
			numeric_vector tmp(m_systemVector.block(i));
			tmp.block_write(file);
		}
#else
		m_systemVector.block_write(file);
#endif
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while writing the system vector to file: Please make sure you have writing permission. Quick fix: Remove StreamingMatrices from OutputFlags");
	}
}
template void SEDGMinLee<2>::saveMatricesToFiles(const string& directory) const;
template void SEDGMinLee<3>::saveMatricesToFiles(const string& directory) const;

template<size_t dim>
void SEDGMinLee<dim>::loadMatricesFromFiles(const string& directory) {
// read the system matrices from file
	try {
		for (size_t i = 0; i < m_stencil->getQ() - 1; i++) {
			for (size_t j = 0; j < m_stencil->getQ() - 1; j++) {
				// filename
				std::stringstream filename;
				filename << directory << "/checkpoint_system_matrix_" << i
#ifdef WITH_TRILINOS_MPI
						<< "_" << j << "_" << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) << ".dat";
#else
						<< "_" << j << ".dat";
#endif
				std::ifstream file(filename.str().c_str());
#ifndef WITH_TRILINOS
				m_systemMatrix.block(i, j).block_read(file);
#else
				// TODO use trilinos functions for read and write. This here is really bad
#endif
			}
		}

	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while reading the system matrices from file: Please switch off the restart option to start the simulation from the beginning.");
	}

// Read the system vector
	try {
		// filename
		std::stringstream filename;
		filename << directory << "/checkpoint_system_vector.dat";
		std::ifstream file(filename.str().c_str());
#ifdef WITH_TRILINOS
		for (size_t i = 0; i < m_stencil->getQ() - 1; i++) {
			numeric_vector tmp(m_systemVector.block(i));
			tmp.block_read(file);
			m_systemVector.block(i) = tmp;
		}
#else

		m_systemVector.block_read(file);
#endif
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		throw AdvectionSolverException(
				"An error occurred while reading the systemVector from file: Please switch off the restart option to start the simulation from the beginning.");
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
	outfile << m_stencil->getD() << endl;
//write Q
	outfile << m_stencil->getQ() << endl;
//write magic number of cell geometry
	outfile << calcMagicNumber() << endl;
//write dqScaling1
	outfile << m_stencil->getDirection(1)(0) << endl;
//write dqScaling2
	outfile << m_stencil->getDirection(1)(1) << endl;
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
	if (tmp != m_stencil->getD()) {
		message = "Dimension not equal.";
		return false;
	}
// Q
	infile >> tmp;
	if (tmp != m_stencil->getQ()) {
		message = "Number of particle velocities not equal.";
		return false;
	}
// magic number of cell geometry
	double dtmp;
	infile >> dtmp;
	if (fabs(dtmp - calcMagicNumber()) > 1e-1) {
		message = "Mesh (or at least its magic number) not equal.";
		return false;
	}
// dqScaling1
	infile >> dtmp;
	if (fabs(dtmp - m_stencil->getDirection(1)(0)) > 1e-5) {
		message = "Scaling of Stencil (1st coordinate) not equal.";
		return false;
	}
// dqScaling2
	infile >> dtmp;
	if (fabs(dtmp - m_stencil->getDirection(1)(1)) > 1e-5) {
		message = "Scaling of Stencil (2nd) not equal.";
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
void SEDGMinLee<dim>::loadCheckpoint(const string& directory) {
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
void SEDGMinLee<dim>::saveCheckpoint(const string& directory) const {
	writeStatus(directory);
	saveMatricesToFiles(directory);
}
template void SEDGMinLee<2>::saveCheckpoint(const string& directory) const;
template void SEDGMinLee<3>::saveCheckpoint(const string& directory) const;

} /* namespace natrium */
