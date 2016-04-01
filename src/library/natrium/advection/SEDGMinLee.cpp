/**
 * @file SEDGMinLee.cpp
 * @short Advection scheme proposed by Min and Lee (2011)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SEDGMinLee.h"

#include <fstream>

#ifdef WITH_TRILINOS
#include "deal.II/lac/trilinos_sparsity_pattern.h"
#else
#include "deal.II/lac/dynamic_sparsity_pattern.h"
#endif

#include "deal.II/dofs/dof_renumbering.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/lac/matrix_iterator.h"
#include "deal.II/lac/sparsity_tools.h"
#include "deal.II/base/utilities.h"
#include "deal.II/lac/constraint_matrix.h"

#include "../problemdescription/PeriodicBoundary.h"
#include "../problemdescription/LinearBoundary.h"

#include "../stencils/Stencil.h"

#include "../utilities/DealiiExtensions.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
SEDGMinLee<dim>::SEDGMinLee(boost::shared_ptr<Mesh<dim> > triangulation,
		boost::shared_ptr<BoundaryCollection<dim> > boundaries,
		size_t orderOfFiniteElement, boost::shared_ptr<Stencil> Stencil,
		bool useCentralFlux) :
		m_mesh(triangulation), m_boundaries(boundaries), m_mapping(
				orderOfFiniteElement), m_stencil(Stencil), m_orderOfFiniteElement(
				orderOfFiniteElement), m_useCentralFlux(useCentralFlux) {
	// assertions
	assert(orderOfFiniteElement >= 1);
	assert(Stencil->getD() == dim);

	// make dof handler
	m_quadrature = boost::make_shared<QGaussLobatto<dim> >(
			orderOfFiniteElement + 1);
	m_faceQuadrature = boost::make_shared<QGaussLobatto<dim - 1> >(
			orderOfFiniteElement + 1);
	m_fe = boost::make_shared<FE_DGQArbitraryNodes<dim> >(
			QGaussLobatto<1>(orderOfFiniteElement + 1));
	m_doFHandler = boost::make_shared<DoFHandler<dim> >(*triangulation);

} /* SEDGMinLee<dim>::SEDGMinLee */

template<size_t dim>
void SEDGMinLee<dim>::setupDoFs() {

	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);

	updateSparsityPattern();

	// define relation between dofs and quadrature nodes
	m_facedof_to_q_index = map_facedofs_to_q_index();
	m_celldof_to_q_index = map_celldofs_to_q_index();
	m_q_index_to_facedof = map_q_index_to_facedofs();

	// set size for the system vector	
	m_systemVector.reinit(m_stencil->getQ() - 1);
	for (size_t i = 0; i < m_stencil->getQ() - 1; i++) {
		m_systemVector.block(i).reinit(m_locallyOwnedDofs, MPI_COMM_WORLD);
		m_systemVector.collect_sizes();
	}

}

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
			| update_quadrature_points | update_JxW_values
			| update_normal_vectors;
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
	std::vector<dealii::types::global_dof_index> localDoFIndices(dofs_per_cell);

///////////////
// MAIN LOOP //
///////////////
	typename DoFHandler<dim>::active_cell_iterator cell =
			m_doFHandler->begin_active(), endc = m_doFHandler->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// calculate the fe values for the cell
			feCellValues.reinit(cell);

			// get global degrees of freedom
			cell->get_dof_indices(localDoFIndices);

			// assemble local cell matrices
			assembleLocalMassMatrix(feCellValues, dofs_per_cell,
					localMassMatrix);
			assembleLocalDerivativeMatrices(feCellValues, dofs_per_cell,
					localDerivativeMatrices);

			// invert local mass matrix
			for (size_t i = 0; i < dofs_per_cell; i++) {
				inverseLocalMassMatrix.at(i) = 1. / localMassMatrix.at(i);
			}

			// assemble faces and put together
			for (size_t alpha = 1; alpha < m_stencil->getQ(); alpha++) {
// calculate local diagonal block (cell) matrix -D
				calculateAndDistributeLocalStiffnessMatrix(alpha,
						localDerivativeMatrices, localSystemMatrix,
						inverseLocalMassMatrix, localDoFIndices, dofs_per_cell);
// calculate face contributions  R
				assembleAndDistributeLocalFaceMatrices(alpha, cell,
						feFaceValues, feSubfaceValues, feNeighborFaceValues,
						inverseLocalMassMatrix);
			}
		}
	}
	m_systemVector.compress(dealii::VectorOperation::add);
	m_systemMatrix.compress(dealii::VectorOperation::add);

//#endif

} /* reassemble */

template<size_t dim>
void SEDGMinLee<dim>::updateSparsityPattern() {
	// TODO only update sparsity pattern for changed cells
	///////////////////////////////////////////////////////
	// Setup sparsity pattern (completely manually):
	///////////////////////////////////////////////////////
	// allocate sizes
	size_t n_blocks = m_stencil->getQ() - 1;
	const dealii::UpdateFlags faceUpdateFlags = update_values
			| update_quadrature_points;
	dealii::FEFaceValues<dim>* feFaceValues = new FEFaceValues<dim>(m_mapping,
			*m_fe, *m_faceQuadrature, faceUpdateFlags);

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

	//get locally owned and locally relevant dofs
	m_locallyOwnedDofs = m_doFHandler->locally_owned_dofs();
	DoFTools::extract_locally_relevant_dofs(*m_doFHandler,
			m_locallyRelevantDofs);

	TrilinosWrappers::SparsityPattern cSparseDiag(m_locallyOwnedDofs,
			m_locallyOwnedDofs, m_locallyRelevantDofs, MPI_COMM_WORLD);
	TrilinosWrappers::SparsityPattern cSparseOpposite(m_locallyOwnedDofs,
			m_locallyOwnedDofs, m_locallyRelevantDofs, MPI_COMM_WORLD);
	TrilinosWrappers::SparsityPattern cSparseNotOpposite(m_locallyOwnedDofs,
			m_locallyOwnedDofs, m_locallyRelevantDofs, MPI_COMM_WORLD);
	/*DynamicSparsityPattern cSparseDiag(m_locallyRelevantDofs);
	 DynamicSparsityPattern cSparseOpposite(m_locallyRelevantDofs);
	 DynamicSparsityPattern cSparseEmpty(m_locallyRelevantDofs);*/

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
	ConstraintMatrix constraints;
	DealIIExtensions::make_sparser_flux_sparsity_pattern(*m_doFHandler,
			cSparseDiag, constraints, *m_boundaries, feFaceValues, true,
			dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
	delete feFaceValues;
	/*DoFTools::make_flux_sparsity_pattern(*m_doFHandler,
	 cSparse.block(0, 0));*/

	// add entries for non-periodic boundaries
	for (typename BoundaryCollection<dim>::ConstLinearIterator dirichlet_iterator =
			m_boundaries->getLinearBoundaries().begin();
			dirichlet_iterator != m_boundaries->getLinearBoundaries().end();
			dirichlet_iterator++) {
		dirichlet_iterator->second->addToSparsityPattern(cSparseOpposite,
				*m_doFHandler);
		if (BoundaryTools::COUPLE_ALL_DISTRIBUTIONS
				== dirichlet_iterator->second->m_distributionCoupling) {
			dirichlet_iterator->second->addToSparsityPattern(cSparseNotOpposite,
					*m_doFHandler);
		}
	}
	//reinitialize matrices
	//In order to store the sparsity pattern for blocks with same pattern only once: initialize from other block
	cSparseDiag.compress();
	cSparseOpposite.compress();
	cSparseNotOpposite.compress();
	/*SparsityTools::distribute_sparsity_pattern(cSparseDiag,
	 m_doFHandler->n_locally_owned_dofs_per_processor(), MPI_COMM_WORLD,
	 m_locallyRelevantDofs);
	 SparsityTools::distribute_sparsity_pattern(cSparseOpposite,
	 m_doFHandler->n_locally_owned_dofs_per_processor(), MPI_COMM_WORLD,
	 m_locallyRelevantDofs);
	 SparsityTools::distribute_sparsity_pattern(cSparseEmpty,
	 m_doFHandler->n_locally_owned_dofs_per_processor(), MPI_COMM_WORLD,
	 m_locallyRelevantDofs);
	 */
	m_systemMatrix.reinit(n_blocks, n_blocks);
	size_t first_opposite = m_stencil->getIndexOfOppositeDirection(1) - 1;
	size_t some_nonopposite = m_stencil->getIndexOfOppositeDirection(1);
	assert(some_nonopposite <= n_blocks);
	assert(some_nonopposite != first_opposite);
	assert(some_nonopposite != 1);
	m_systemMatrix.block(0, 0).reinit(cSparseDiag);
	m_systemMatrix.block(0, some_nonopposite).reinit(cSparseNotOpposite);
	m_systemMatrix.block(0, first_opposite).reinit(cSparseOpposite);

	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			if ((I == 0) and (J == 0)) {
				continue;
			}
			if ((I == 0) and (J == some_nonopposite)) {
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
				m_systemMatrix.block(I, J).reinit(
						m_systemMatrix.block(0, first_opposite));
				continue;
			} else {
				m_systemMatrix.block(I, J).reinit(
						m_systemMatrix.block(0, some_nonopposite));
			}

		}
	}
	m_systemMatrix.collect_sizes();

}
/* updateSparsityPattern */

template<size_t dim>
void SEDGMinLee<dim>::assembleLocalMassMatrix(
		const dealii::FEValues<dim>& feValues, size_t dofs_per_cell,
		vector<double> &massMatrix) {
// initialize with zeros
	std::fill(massMatrix.begin(), massMatrix.end(), 0.0);

// fill diagonal "matrix"
	for (size_t i = 0; i < dofs_per_cell; i++) {
		size_t q_point = m_celldof_to_q_index.at(i);
		massMatrix.at(i) += feValues.shape_value(i, q_point)
				* feValues.shape_value(i, q_point) * feValues.JxW(q_point);
	}

} /*assembleLocalMassMatrix*/

template<size_t dim>
void SEDGMinLee<dim>::assembleLocalDerivativeMatrices(
		const dealii::FEValues<dim>& feValues, size_t dofs_per_cell,
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

template<size_t dim>
void SEDGMinLee<dim>::assembleAndDistributeLocalFaceMatrices(size_t alpha,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		dealii::FEFaceValues<dim>& feFaceValues,
		dealii::FESubfaceValues<dim>& feSubfaceValues,
		dealii::FEFaceValues<dim>& feNeighborFaceValues,
		const vector<double>& inverseLocalMassMatrix) {

// loop over all faces
	for (size_t j = 0; j < dealii::GeometryInfo<dim>::faces_per_cell; j++) {
		//Faces at boundary
		if (cell->face(j)->at_boundary()) {
			size_t boundaryIndicator = cell->face(j)->boundary_id();
			if (m_boundaries->isPeriodic(boundaryIndicator)) {
				// Apply periodic boundaries
				const boost::shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
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
				if ((m_boundaries->getBoundary(boundaryIndicator)->isLinear())) {
					const boost::shared_ptr<LinearBoundary<dim> >& LinearBoundary =
							m_boundaries->getLinearBoundary(boundaryIndicator);
					LinearBoundary->assembleBoundary(alpha, cell, j,
							feFaceValues, *m_stencil,
							m_q_index_to_facedof.at(j), inverseLocalMassMatrix,
							m_systemMatrix, m_systemVector, m_useCentralFlux);
				}
			} /* endif isPeriodic */

		} else /* if is not at boundary */{
// Internal faces
			typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(
					j);
			assembleAndDistributeInternalFace(alpha, cell, j, neighbor,
					cell->neighbor_face_no(j), feFaceValues, feSubfaceValues,
					feNeighborFaceValues, inverseLocalMassMatrix);
		} /* endif is face at boundary*/

	}

} /* assembleLocalFaceMatrix */

template<> void SEDGMinLee<2>::calculateAndDistributeLocalStiffnessMatrix(
		size_t alpha,
		const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix,
		const vector<double>& inverseLocalMassMatrix,
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
			block.add(globalDoFs[j], globalDoFs[k],
					systemMatrix(j, k) * inverseLocalMassMatrix.at(j));
		}
}
template<> void SEDGMinLee<3>::calculateAndDistributeLocalStiffnessMatrix(
		size_t alpha,
		const vector<dealii::FullMatrix<double> >& derivativeMatrices,
		dealii::FullMatrix<double> &systemMatrix,
		const vector<double>& inverseLocalMassMatrix,
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
			block.add(globalDoFs[j], globalDoFs[k],
					systemMatrix(j, k) * inverseLocalMassMatrix.at(j));
}

template<size_t dim>
void SEDGMinLee<dim>::assembleAndDistributeInternalFace(size_t alpha,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber,
		typename dealii::DoFHandler<dim>::cell_iterator& neighborCell,
		size_t neighborFaceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		dealii::FESubfaceValues<dim>& feSubfaceValues,
		dealii::FEFaceValues<dim>& feNeighborFaceValues,
		const vector<double>& inverseLocalMassMatrix) {
// get the required FE Values for the local cell
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<Tensor<1, dim> > &normals =
			feFaceValues.get_all_normal_vectors();

	if (4 == alpha) {

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
		assert(feFaceValues.shape_value(thisDoF, q) > 0);
		size_t neighborDoF = m_q_index_to_facedof.at(neighborFaceNumber).at(q);
		assert(feNeighborFaceValues.shape_value(neighborDoF, q) > 0);

		double cell_entry = 0.0;
		double neighbor_entry = 0.0;

		// calculate matrix entries
		double prefactor = JxW.at(q);
		double exn = 0.0;

		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {		// TODO efficient multiplication
			exn += normals.at(q)[i] * m_stencil->getDirection(alpha)(i);
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

template<size_t dim>
void SEDGMinLee<dim>::stream() {
}

template class SEDGMinLee<2> ;
template class SEDGMinLee<3> ;

} /* namespace natrium */
