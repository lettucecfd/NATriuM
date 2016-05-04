/**
 * @file SemiLagrangian.cpp
 * @short Semi-Lagrangian advection scheme
 * @date 29.04.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SemiLagrangian.h"

#include <fstream>

#include "deal.II/lac/trilinos_sparsity_pattern.h"
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
#include "../utilities/CFDSolverUtilities.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
SemiLagrangian<dim>::SemiLagrangian(boost::shared_ptr<Mesh<dim> > triangulation,
		boost::shared_ptr<BoundaryCollection<dim> > boundaries,
		size_t orderOfFiniteElement, boost::shared_ptr<Stencil> Stencil,
		double delta_t) :
		m_mesh(triangulation), m_boundaries(boundaries), m_mapping(
				orderOfFiniteElement), m_stencil(Stencil), m_orderOfFiniteElement(
				orderOfFiniteElement), m_deltaT(delta_t) {
	// assertions
	assert(orderOfFiniteElement >= 1);
	assert(Stencil->getD() == dim);

	// TODO FEM rather than DG
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
void SemiLagrangian<dim>::setupDoFs() {

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
void SemiLagrangian<dim>::reassemble() {
// TODO: if Mesh changed: reinit dof-handler and sparsity pattern in some way

	// make sure that sparsity structure is not empty
	assert(m_systemMatrix.n() != 0);
	assert(m_systemMatrix.m() != 0);

/////////////////////////////////
// Initialize Finite Element ////
/////////////////////////////////
// Define update flags (which values have to be known at each cell, face, neighbor face)
	const dealii::UpdateFlags cellUpdateFlags = update_values
			| update_quadrature_points | update_JxW_values;
	const dealii::UpdateFlags faceUpdateFlags = update_quadrature_points
			| update_JxW_values | update_normal_vectors;
// Finite Element
	dealii::FEValues<dim> feCellValues(m_mapping, *m_fe, *m_quadrature,
			cellUpdateFlags);
	dealii::FEFaceValues<dim> feFaceValues(m_mapping, *m_fe, *m_faceQuadrature,
			faceUpdateFlags);
	dealii::FESubfaceValues<dim> feSubfaceValues(m_mapping, *m_fe,
			*m_faceQuadrature, faceUpdateFlags);

	const size_t dofs_per_cell = m_fe->dofs_per_cell;

// Initialize

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

		}
	}
	m_systemVector.compress(dealii::VectorOperation::add);
	m_systemMatrix.compress(dealii::VectorOperation::add);

//#endif

} /* reassemble */

template<size_t dim>
void SemiLagrangian<dim>::updateSparsityPattern() {
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
std::map<size_t, size_t> SemiLagrangian<dim>::map_celldofs_to_q_index() const {
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
vector<std::map<size_t, size_t> > SemiLagrangian<dim>::map_facedofs_to_q_index() const {
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
vector<std::map<size_t, size_t> > SemiLagrangian<dim>::map_q_index_to_facedofs() const {
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
void SemiLagrangian<dim>::stream() {
}

template<size_t dim>
void SemiLagrangian<dim>::fillSparseObject(bool sparsity_pattern) {

	if (not sparsity_pattern) {
		// make sure that sparsity structure is not empty
		assert(m_systemMatrix.n() != 0);
		assert(m_systemMatrix.m() != 0);

	}

	/////////////////////////////////
	// Initialize Finite Element ////
	/////////////////////////////////
	// Define update flags (which values have to be known at each cell, face, neighbor face)
	const dealii::UpdateFlags cell_update_flags = update_values
			| update_quadrature_points | update_JxW_values;
	const dealii::UpdateFlags face_update_flags = update_quadrature_points
			| update_normal_vectors;
	// Finite Element
	dealii::FEValues<dim> fe_cell_values(m_mapping, *m_fe, *m_quadrature,
			cell_update_flags);
	dealii::FEFaceValues<dim> fe_face_values(m_mapping, *m_fe,
			*m_faceQuadrature, face_update_flags);

	// Initialize
	std::map<typename DoFHandler<dim>::active_cell_iterator, XList> cell_map;
	const size_t dofs_per_cell = m_fe->dofs_per_cell;
	const std::vector<Point<dim> > & unit_support_points =
			m_fe->get_unit_support_points();
	std::vector<Tensor<1, dim> > minus_dtealpha;
	for (size_t i = 0; i < m_stencil->getQ(); i++) {
		minus_dtealpha.push_back(
				vectorToTensor(m_stencil->getDirection(i)) * (-m_deltaT));
	}
	std::vector<dealii::types::global_dof_index> local_dof_indices(
			dofs_per_cell);
	std::vector<dealii::types::global_dof_index> neighbor_dof_indices(
			dofs_per_cell);
	size_t max_n_shells = (size_t) (m_stencil->getMaxParticleVelocityMagnitude()
			* m_deltaT
			/ CFDSolverUtilities::getMinimumVertexDistance<dim>(*m_mesh)) + 1;
	if ((max_n_shells > 1)
			and (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1)) {
		LOG(WARNING)
				<< "The global CFL number (wrt. cells) is > 1. That is not a problem, in principle. "
						"However, depending on the mesh and its distribution, the simulation may crash, as the "
						"semi-Lagrangian paths may lead into cells that are not even ghost cells. "
						"When this happens, the algorithm is not longer defined and you will get an error message "
						"before the simulation crashes.";
	}

	///////////////
	// MAIN LOOP //
	///////////////
	typename DoFHandler<dim>::active_cell_iterator cell =
			m_doFHandler->begin_active(), endc = m_doFHandler->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// calculate the fe values for the cell
			fe_cell_values.reinit(cell);

			//
			cell_map.clear();

			// get a list of adjacent cells (including all cells that have common vertices with 'cell')
			Neighborhood neighborhood;
			getNeighborhood(cell, neighborhood, max_n_shells);

			// get global degrees of freedom
			cell->get_dof_indices(local_dof_indices);

			// for all points in cell
			for (size_t i = 0; i < dofs_per_cell; i++) {
				// get a point x
				dealii::Point<dim> x_i = m_mapping.transform_unit_to_real_cell(
						cell, unit_support_points.at(i));
				// for all directions
				for (size_t alpha = 1; alpha < m_stencil->getQ(); alpha++) {
					// calculate x^(t-delta_t)
					dealii::Point<dim> x_previous = x_i
							+ minus_dtealpha.at(alpha);
					bool found = false;
					DoFInfo info_i(local_dof_indices.at(i), alpha, alpha,
							x_previous);

					// look for x_previous in neighborhood
					for (size_t j = 0; j < neighborhood.size(); j++) {
						if (neighborhood.at(j)->point_inside(x_previous)) {
							found = true;
							typename std::map<
									typename DoFHandler<dim>::active_cell_iterator, XList>::iterator it =
									cell_map.find(neighborhood.at(j));
							if (it == cell_map.end()) {
								XList new_xlist;
								new_xlist.push_back(info_i);
								cell_map.insert(
										std::pair<
												typename DoFHandler<dim>::active_cell_iterator,
												XList>(neighborhood.at(j),
												new_xlist));
							} else {
								it->second.push_back(info_i);
							}
							break;
						}
					}
				}
			}
		}

	}
	m_systemVector.compress(dealii::VectorOperation::add);
	m_systemMatrix.compress(dealii::VectorOperation::add);

//#endif
} /* fillSparseObject */

template class SemiLagrangian<2> ;
template class SemiLagrangian<3> ;

} /* namespace natrium */
