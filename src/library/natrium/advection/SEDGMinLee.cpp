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
#include <deal.II/lac/affine_constraints.h>

#include "../boundaries/PeriodicBoundary.h"
#include "../boundaries/Boundary.h"
#include "../utilities/NATriuMException.h"

#include "../stencils/Stencil.h"

#include "../utilities/DealiiExtensions.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
SEDGMinLee<dim>::SEDGMinLee(ProblemDescription<dim>& problem,
		size_t orderOfFiniteElement, QuadratureName quad_name,
		SupportPointsName points_name, boost::shared_ptr<Stencil> stencil, bool use_central_flux,
		double delta_t) :
		AdvectionOperator<dim>(problem, orderOfFiniteElement, quad_name,
				points_name, stencil, delta_t, true), m_useCentralFlux(
				use_central_flux) {
	if (use_central_flux){
		LOG(WARNING) << "The use of central fluxes is highly discouraged. They "
				"have not proven too good and have not been used or tested for a long time." << endl;
	}
	if ((quad_name != QGAUSS_LOBATTO) or (points_name != GAUSS_LOBATTO_POINTS) ){
		throw AdvectionSolverException("SEDGMinLee can only be used with Gauss-Lobatto points"
				"and Gauss-Lobatto quadrature. Everything else would probably be way too costly and is thus not"
				"implemented in NATriuM.");
	}

} /* SEDGMinLee<dim>::SEDGMinLee */

template<size_t dim>
SEDGMinLee<dim>::SEDGMinLee(ProblemDescription<dim>& problem,
		size_t orderOfFiniteElement, boost::shared_ptr<Stencil> stencil,
		double delta_t) :
		SEDGMinLee(problem, orderOfFiniteElement, QGAUSS_LOBATTO,
				GAUSS_LOBATTO_POINTS, stencil, delta_t, false) {

}

template<size_t dim>
void SEDGMinLee<dim>::setupDoFs() {

	// distribute degrees of freedom over mesh
	Base::distributeDoFs();

	updateSparsityPattern();

	// define relation between dofs and quadrature nodes
	m_facedof_to_q_index = map_facedofs_to_q_index();
	m_celldof_to_q_index = map_celldofs_to_q_index();
	m_q_index_to_facedof = map_q_index_to_facedofs();

	// set size for the system vector	
	m_systemVector.reinit(Base::m_stencil->getQ() - 1);
	for (size_t i = 0; i < Base::m_stencil->getQ() - 1; i++) {
		m_systemVector.block(i).reinit(Base::getLocallyOwnedDofs(),
		MPI_COMM_WORLD);
		m_systemVector.collect_sizes();
	}

}

template<size_t dim>
void SEDGMinLee<dim>::reassemble() {
// TODO: if Mesh changed: reinit dof-handler and sparsity pattern in some way

// make sure that sparsity structure is not empty
	distributed_sparse_block_matrix& system_matrix = Base::m_systemMatrix;
	assert(system_matrix.n() != 0);
	assert(system_matrix.m() != 0);
/////////////////////////////////
// Initialize Finite Element ////
/////////////////////////////////
// Define update flags (which values have to be known at each cell, face, neighbor face)
	const dealii::UpdateFlags cell_flags = update_values | update_gradients
			| update_quadrature_points | update_JxW_values
			| update_inverse_jacobians;
	const dealii::UpdateFlags face_flags = update_values
			| update_quadrature_points | update_JxW_values
			| update_normal_vectors;
	const dealii::UpdateFlags neighbor_face_flags = update_values
			| update_JxW_values | update_normal_vectors;
// Finite Element
	boost::shared_ptr<dealii::FEValues<dim> > fe_cell_ptr = Base::getFEValues(
			cell_flags);
	dealii::FEValues<dim>& fe_cell_values = *fe_cell_ptr;
	boost::shared_ptr<dealii::FEFaceValues<dim> > fe_face_ptr =
			Base::getFEFaceValues(face_flags);
	dealii::FEFaceValues<dim> & fe_face_values = *fe_face_ptr;
	boost::shared_ptr<dealii::FEFaceValues<dim> > fe_neighbor_face_ptr =
			Base::getFEFaceValues(neighbor_face_flags);
	dealii::FEFaceValues<dim> & fe_neighbor_face_values = *fe_neighbor_face_ptr;
	dealii::FESubfaceValues<dim> fe_subface_values(*Base::m_mapping,
			*Base::m_fe, *Base::m_faceQuadrature, face_flags);
	const size_t dofs_per_cell = Base::m_fe->dofs_per_cell;
	const size_t Q = Base::getStencil()->getQ();

// Initialize matrices
	vector<double> local_mass_matrix(dofs_per_cell);
	vector<double> inverse_local_mass_matrix(dofs_per_cell);
	vector<dealii::FullMatrix<double> > local_stiffness_matrix;
	for (size_t i = 0; i < dim; i++) {
		dealii::FullMatrix<double> D_i(dofs_per_cell, dofs_per_cell);
		local_stiffness_matrix.push_back(D_i);
	}
	dealii::FullMatrix<double> local_face_matrix(dofs_per_cell, dofs_per_cell);
	dealii::FullMatrix<double> local_system_matrix(dofs_per_cell,
			dofs_per_cell);
	std::vector<dealii::types::global_dof_index> local_doF_indices(
			dofs_per_cell);

///////////////
// MAIN LOOP //
///////////////
	typename DoFHandler<dim>::active_cell_iterator cell =
			Base::m_doFHandler->begin_active(), endc =
			Base::m_doFHandler->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// calculate the fe values for the cell
			fe_cell_values.reinit(cell);

			// get global degrees of freedom
			cell->get_dof_indices(local_doF_indices);

			// assemble local cell matrices
			assembleLocalMassMatrix(fe_cell_values, dofs_per_cell,
					local_mass_matrix);
			assembleLocalDerivativeMatrices(fe_cell_values, dofs_per_cell,
					local_stiffness_matrix);

			// invert local mass matrix
			for (size_t i = 0; i < dofs_per_cell; i++) {
				inverse_local_mass_matrix.at(i) = 1. / local_mass_matrix.at(i);
			}

			// assemble faces and put together
			for (size_t alpha = 1; alpha < Q; alpha++) {
// calculate local diagonal block (cell) matrix -D
				calculateAndDistributeLocalStiffnessMatrix(alpha,
						local_stiffness_matrix, local_system_matrix,
						inverse_local_mass_matrix, local_doF_indices,
						dofs_per_cell);
// calculate face contributions  R
				assembleAndDistributeLocalFaceMatrices(alpha, cell,
						fe_face_values, fe_subface_values,
						fe_neighbor_face_values, inverse_local_mass_matrix);
			}
		}
	}
	m_systemVector.compress(dealii::VectorOperation::add);
	Base::m_systemMatrix.compress(dealii::VectorOperation::add);

//#endif

} /* reassemble */

template<size_t dim>
void SEDGMinLee<dim>::updateSparsityPattern() {
	// TODO only update sparsity pattern for changed cells
	///////////////////////////////////////////////////////
	// Setup sparsity pattern (completely manually):
	///////////////////////////////////////////////////////
	// allocate sizes
	size_t n_blocks = Base::m_stencil->getQ() - 1;
	const dealii::UpdateFlags face_flags = update_values
			| update_quadrature_points;
	dealii::FEFaceValues<dim>* fe_face_values = new FEFaceValues<dim>(
			*Base::m_mapping, *Base::m_fe, *Base::m_faceQuadrature, face_flags);

	// create cell maps for periodic boundary
	for (typename BoundaryCollection<dim>::ConstPeriodicIterator periodic =
			Base::getBoundaries()->getPeriodicBoundaries().begin();
			periodic != Base::getBoundaries()->getPeriodicBoundaries().end();
			periodic++) {
		// Periodic boundaries have two boundary indicators; (and are stored twice in the map)
		// skip double execution of addToSparsityPattern
		if (periodic->first == periodic->second->getBoundaryIndicator1()) {
			periodic->second->createCellMap(*Base::m_doFHandler);
		}
	}

	TrilinosWrappers::SparsityPattern cSparseDiag(Base::getLocallyOwnedDofs(),
			Base::getLocallyOwnedDofs(), Base::getLocallyRelevantDofs(),
			MPI_COMM_WORLD);
	TrilinosWrappers::SparsityPattern cSparseOpposite(
			Base::getLocallyOwnedDofs(), Base::getLocallyOwnedDofs(),
			Base::getLocallyRelevantDofs(), MPI_COMM_WORLD);
	TrilinosWrappers::SparsityPattern cSparseNotOpposite(
			Base::getLocallyOwnedDofs(), Base::getLocallyOwnedDofs(),
			Base::getLocallyRelevantDofs(), MPI_COMM_WORLD);
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
	AffineConstraints<double> constraints;
	DealIIExtensions::make_sparser_flux_sparsity_pattern(*Base::m_doFHandler,
			cSparseDiag, constraints, *Base::getBoundaries(), fe_face_values,
			true, dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
	delete fe_face_values;
	/*DoFTools::make_flux_sparsity_pattern(*m_doFHandler,
	 cSparse.block(0, 0));*/

	// add entries for non-periodic boundaries
	for (typename BoundaryCollection<dim>::ConstIterator boundary_it =
			Base::getBoundaries()->getBoundaries().begin();
			boundary_it != Base::getBoundaries()->getBoundaries().end();
			boundary_it++) {
		if (boundary_it->second->isPeriodic())
			continue;
		if (not boundary_it->second->isDGSupported()){
			throw NotImplementedException("SEDG so far only supports boundaries that"
					" are linear in f. A boundary from your ProblemDescriptions"
					" had isDGSupported() == false.");
		}
		boundary_it->second->addToSparsityPattern(cSparseOpposite,
				*Base::m_doFHandler);
		if (BoundaryTools::COUPLE_ALL_DISTRIBUTIONS
				== boundary_it->second->getDistributionCoupling()) {
			boundary_it->second->addToSparsityPattern(cSparseNotOpposite,
					*Base::m_doFHandler);
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
	Base::m_systemMatrix.reinit(n_blocks, n_blocks);
	size_t first_opposite = Base::m_stencil->getIndexOfOppositeDirection(1) - 1;
	size_t some_nonopposite = Base::m_stencil->getIndexOfOppositeDirection(1);
	assert(some_nonopposite <= n_blocks);
	assert(some_nonopposite != first_opposite);
	assert(some_nonopposite != 1);
	Base::m_systemMatrix.block(0, 0).reinit(cSparseDiag);
	Base::m_systemMatrix.block(0, some_nonopposite).reinit(cSparseNotOpposite);
	Base::m_systemMatrix.block(0, first_opposite).reinit(cSparseOpposite);

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
				Base::m_systemMatrix.block(I, J).reinit(
						Base::m_systemMatrix.block(0, 0));
				continue;
			}
			if (I == Base::m_stencil->getIndexOfOppositeDirection(J + 1) - 1) {
				Base::m_systemMatrix.block(I, J).reinit(
						Base::m_systemMatrix.block(0, first_opposite));
				continue;
			} else {
				Base::m_systemMatrix.block(I, J).reinit(
						Base::m_systemMatrix.block(0, some_nonopposite));
			}

		}
	}
	Base::m_systemMatrix.collect_sizes();

}
/* updateSparsityPattern */

template<size_t dim>
void SEDGMinLee<dim>::assembleLocalMassMatrix(
		const dealii::FEValues<dim>& fe_values, size_t dofs_per_cell,
		vector<double> &mass_matrix) {
// initialize with zeros
	std::fill(mass_matrix.begin(), mass_matrix.end(), 0.0);

// fill diagonal "matrix"
	for (size_t i = 0; i < dofs_per_cell; i++) {
		size_t q_point = m_celldof_to_q_index.at(i);
		mass_matrix.at(i) += fe_values.shape_value(i, q_point)
				* fe_values.shape_value(i, q_point) * fe_values.JxW(q_point);
	}

} /*assembleLocalMassMatrix*/

template<size_t dim>
void SEDGMinLee<dim>::assembleLocalDerivativeMatrices(
		const dealii::FEValues<dim>& fe_values, size_t dofs_per_cell,
		vector<dealii::FullMatrix<double> >&derivative_matrix) const {
	for (size_t i = 0; i < dim; i++) {
		derivative_matrix.at(i) = 0;
	}
	for (size_t i = 0; i < dofs_per_cell; i++) {
		for (size_t j = 0; j < dofs_per_cell; j++) {
			// the shape value is zero for all q, except i<->q
			size_t q_point = m_celldof_to_q_index.at(i);
			Tensor<1, dim> integrandAtQ;
			integrandAtQ = fe_values.shape_grad(j, q_point);
			integrandAtQ *= (fe_values.shape_value(i, q_point)
					* fe_values.JxW(q_point));
			for (size_t k = 0; k < dim; k++) {
				derivative_matrix.at(k)(i, j) += integrandAtQ[k];
			}

		}
	}

} /* assembleLocalDerivativeMatrix */

template<size_t dim>
void SEDGMinLee<dim>::assembleAndDistributeLocalFaceMatrices(size_t alpha,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		dealii::FEFaceValues<dim>& fe_face_values,
		dealii::FESubfaceValues<dim>& fe_subface_values,
		dealii::FEFaceValues<dim>& fe_neighbor_face_values,
		const vector<double>& inverse_local_mass_matrix) {

// loop over all faces
	for (size_t j = 0; j < dealii::GeometryInfo<dim>::faces_per_cell; j++) {
		//Faces at boundary
		if (cell->face(j)->at_boundary()) {
			size_t boundary_id = cell->face(j)->boundary_id();
			if (Base::getBoundaries()->isPeriodic(boundary_id)) {
				// Apply periodic boundaries
				const boost::shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
						Base::getBoundaries()->getPeriodicBoundary(boundary_id);
				assert(periodicBoundary->isFaceInBoundary(cell, j));
				typename dealii::DoFHandler<dim>::cell_iterator neighborCell;
				size_t opposite_face =
						periodicBoundary->getOppositeCellAtPeriodicBoundary(
								cell, neighborCell);
				assembleAndDistributeInternalFace(alpha, cell, j, neighborCell,
						opposite_face, fe_face_values, fe_subface_values,
						fe_neighbor_face_values, inverse_local_mass_matrix);
			} else /* if is not periodic */ {
				// Apply other boundaries
				if ((Base::getBoundaries()->getBoundary(boundary_id)->isDGSupported())) {
					Base::getBoundaries()->getBoundary(boundary_id)->assembleBoundary(alpha, cell, j,
							fe_face_values, *Base::m_stencil,
							m_q_index_to_facedof.at(j),
							inverse_local_mass_matrix, Base::m_systemMatrix,
							m_systemVector, m_useCentralFlux);
				} else {
					throw NotImplementedException("The SEDG solver so far only supports boundaries"
							" that are linear in the distribution functions, but your ProblemDescription "
							"includes a boundary that is neither periodic, nor isDGSupported() = true.");
				}
			} /* endif isPeriodic */

		} else /* if is not at boundary */{
// Internal faces
			typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(
					j);
			assembleAndDistributeInternalFace(alpha, cell, j, neighbor,
					cell->neighbor_face_no(j), fe_face_values,
					fe_subface_values, fe_neighbor_face_values,
					inverse_local_mass_matrix);
		} /* endif is face at boundary*/

	}

} /* assembleLocalFaceMatrix */

template<> void SEDGMinLee<2>::calculateAndDistributeLocalStiffnessMatrix(
		size_t alpha,
		const vector<dealii::FullMatrix<double> >& derivative_matrices,
		dealii::FullMatrix<double> &system_matrix,
		const vector<double>& inverse_local_mass_matrix,
		const std::vector<dealii::types::global_dof_index>& global_dofs,
		size_t dofs_per_cell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	system_matrix = derivative_matrices.at(0);
	system_matrix *= (-(Base::m_stencil->getDirection(alpha)[0]));
	system_matrix.add(-(Base::m_stencil->getDirection(alpha)[1]),
			derivative_matrices.at(1));
// distribute to global system matrix
	// avoid to call block() too often
	distributed_sparse_matrix& block = Base::m_systemMatrix.block(alpha - 1,
			alpha - 1);
	for (unsigned int j = 0; j < dofs_per_cell; j++)
		for (unsigned int k = 0; k < dofs_per_cell; k++) {
			block.add(global_dofs[j], global_dofs[k],
					system_matrix(j, k) * inverse_local_mass_matrix.at(j));
		}
}
template<> void SEDGMinLee<3>::calculateAndDistributeLocalStiffnessMatrix(
		size_t alpha,
		const vector<dealii::FullMatrix<double> >& derivative_matrices,
		dealii::FullMatrix<double> &system_matrix,
		const vector<double>& inverse_local_mass_matrix,
		const std::vector<dealii::types::global_dof_index>& global_dofs,
		size_t dofs_per_cell) {
// TODO efficient implementation (testing if e_ix, e_iy = 0, -1 or 1)
// calculate -D = -(e_x * D_x  +  e_y * D_y)
	system_matrix = derivative_matrices.at(0);
	system_matrix *= (-Base::m_stencil->getDirection(alpha)[0]);
	system_matrix.add(-Base::m_stencil->getDirection(alpha)[1],
			derivative_matrices.at(1), -Base::m_stencil->getDirection(alpha)[2],
			derivative_matrices.at(2));
// distribute to global system matrix
	distributed_sparse_matrix& block = Base::m_systemMatrix.block(alpha - 1,
			alpha - 1);
	for (unsigned int j = 0; j < dofs_per_cell; j++)
		for (unsigned int k = 0; k < dofs_per_cell; k++)
			block.add(global_dofs[j], global_dofs[k],
					system_matrix(j, k) * inverse_local_mass_matrix.at(j));
}

template<size_t dim>
void SEDGMinLee<dim>::assembleAndDistributeInternalFace(size_t alpha,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t face_number,
		typename dealii::DoFHandler<dim>::cell_iterator& neighbor_cell,
		size_t neighbor_face_number, dealii::FEFaceValues<dim>& fe_face_values,
		dealii::FESubfaceValues<dim>&,
		dealii::FEFaceValues<dim>& fe_neighbor_face_values,
		const vector<double>& inverse_local_mass_matrix) {
// get the required FE Values for the local cell
	fe_face_values.reinit(cell, face_number);
	const vector<double> &JxW = fe_face_values.get_JxW_values();
	const vector<Tensor<1, dim> > &normals =
			fe_face_values.get_normal_vectors();

	if (4 == alpha) {

	}
// get the required dofs of the neighbor cell
//	typename DoFHandler<dim>::face_iterator neighborFace = neighborCell->face(
//			neighborFaceNumber);
	fe_neighbor_face_values.reinit(neighbor_cell, neighbor_face_number);

	vector<dealii::types::global_dof_index> localDoFIndices(
			fe_face_values.get_fe().dofs_per_cell);
	vector<dealii::types::global_dof_index> neighborDoFIndices(
			fe_face_values.get_fe().dofs_per_cell);
	cell->get_dof_indices(localDoFIndices);
	neighbor_cell->get_dof_indices(neighborDoFIndices);

	// avoid many calls to block()
	distributed_sparse_matrix& block = Base::m_systemMatrix.block(alpha - 1,
			alpha - 1);

// loop over all quadrature points at the face
	for (size_t q = 0; q < fe_face_values.n_quadrature_points; q++) {
		size_t thisDoF = m_q_index_to_facedof.at(face_number).at(q);
		assert(fe_face_values.shape_value(thisDoF, q) > 0);
		size_t neighborDoF = m_q_index_to_facedof.at(neighbor_face_number).at(
				q);
		assert(fe_neighbor_face_values.shape_value(neighborDoF, q) > 0);

		double cell_entry = 0.0;
		double neighbor_entry = 0.0;

		// calculate matrix entries
		double prefactor = JxW.at(q);
		double exn = 0.0;

		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {		// TODO efficient multiplication
			exn += normals.at(q)[i] * Base::m_stencil->getDirection(alpha)(i);
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
				cell_entry * inverse_local_mass_matrix.at(thisDoF));
		block.add(localDoFIndices[thisDoF], neighborDoFIndices[neighborDoF],
				neighbor_entry * inverse_local_mass_matrix.at(thisDoF));
	}

// get DoF indices
// TODO cut out construction (allocation); Allocating two vectors in most inner loop is too expensive

// WARNING: Hanging nodes are not implemented, yet
// TODO Implement local refinement

} /* assembleAndDistributeInternalFace */

template<size_t dim>
std::map<size_t, size_t> SEDGMinLee<dim>::map_celldofs_to_q_index() const {
	const dealii::UpdateFlags cell_flags = update_values
			| update_quadrature_points;
// Finite Element
	dealii::FEValues<dim> feCellValues(*Base::m_mapping, *Base::m_fe,
			*Base::m_quadrature, cell_flags);
	const size_t dofs_per_cell = Base::m_fe->dofs_per_cell;
	const size_t n_quadrature_points = Base::m_quadrature->size();
// take first cell
	typename DoFHandler<dim>::active_cell_iterator cell =
			Base::m_doFHandler->begin_active();
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
	const dealii::UpdateFlags face_flags = update_values
			| update_quadrature_points;
	dealii::FEFaceValues<dim> feFaceValues(*Base::m_mapping, *Base::m_fe,
			*Base::m_faceQuadrature, face_flags);

	typename DoFHandler<dim>::active_cell_iterator cell =
			Base::m_doFHandler->begin_active();
// LOOP over all faces
	vector<std::map<size_t, size_t> > result;
	for (size_t f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
		feFaceValues.reinit(cell, f);
		std::map<size_t, size_t> resultForFaceF;
		for (size_t i = 0; i < Base::m_fe->dofs_per_cell; i++) {
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
	dealii::FEFaceValues<dim> feFaceValues(*Base::m_mapping, *Base::m_fe,
			*Base::m_faceQuadrature, faceUpdateFlags);

	typename DoFHandler<dim>::active_cell_iterator cell =
			Base::m_doFHandler->begin_active();
// LOOP over all faces
	vector<std::map<size_t, size_t> > result;
	for (size_t f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
		feFaceValues.reinit(cell, f);
		std::map<size_t, size_t> resultForFaceF;
		for (size_t i = 0; i < Base::m_fe->dofs_per_cell; i++) {
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

template class SEDGMinLee<2> ;
template class SEDGMinLee<3> ;

} /* namespace natrium */
