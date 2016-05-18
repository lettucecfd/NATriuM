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

	fillSparseObject(false); // "false" says that you want to fill the system matrix, not the sparsity pattern

	m_systemVector.compress(dealii::VectorOperation::add);
	m_systemMatrix.compress(dealii::VectorOperation::add);

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

	// create empty sparsity pattern for each block
	for (size_t i = 0; i < n_blocks; i++) {
		std::vector<TrilinosWrappers::SparsityPattern> line;
		line.clear();
		for (size_t j = 0; j < n_blocks; j++) {
			TrilinosWrappers::SparsityPattern empty;
			line.push_back(empty);
		}
		m_sparsityPattern.push_back(line);
	}

	// 	reinit sparsity pattern
	for (size_t i = 0; i < n_blocks; i++)
		for (size_t j = 0; j < n_blocks; j++)
			m_sparsityPattern[i][j].reinit(m_locallyOwnedDofs,
					m_locallyOwnedDofs, m_locallyRelevantDofs, MPI_COMM_WORLD);

	// fill sparsity pattern
	fillSparseObject(true); // "true" means that the object to fill is the sparsity pattern

	// initialize matrix
	m_systemMatrix.reinit(n_blocks, n_blocks);
	for (size_t i = 0; i < n_blocks; i++) {
		for (size_t j = 0; j < n_blocks; j++) {
			m_sparsityPattern[i][j].compress();
			m_systemMatrix.block(i, j).reinit(m_sparsityPattern[i][j]);
		}
	}
	m_systemMatrix.collect_sizes();

	delete feFaceValues;
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
	//dealii::FEValues<dim> fe_normals(m_mappting, *m_fe, *m_quadrature, update_normal_vectors);
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
						"semi-Lagrangian paths may lead into cells that are neither locally owned cells nor ghost cells. "
						"When this happens, the algorithm is not longer well-defined and you will get an error message "
						"before the simulation crashes. The cure: Decrease your time step size (via the CFL number)";
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
									typename DoFHandler<dim>::active_cell_iterator,
									XList>::iterator it = cell_map.find(
									neighborhood.at(j));
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
					} // for j < neighborhood.size()
					if (found)
						continue;

					// calculate x_previous due to boundary conditions
					bool found_across_bound = false;
					int face_hit_counter = -1;
					while (not found_across_bound) {
						face_hit_counter++;
						if (face_hit_counter >= 4) {
							break;
						}
						for (size_t j = 0;
								j < GeometryInfo<dim>::faces_per_cell; j++) {
							if (cell->at_boundary(j)) {
								// check if point is outside
								fe_face_values.reinit(cell, j);
								// assume that all faces are planar (i.e. all normal vectors are equal)
								const dealii::Tensor<1, dim>& n =
										fe_face_values.normal_vector(0);
								const dealii::Point<dim>& p =
										cell->face(j)->vertex(0);
								double n_pminx = n * (p - x_i);
								double n_xpreviousminx = n * (x_previous - x_i);
								if (fabs(n_xpreviousminx) < 1e-20) {
									// Lagrangian path parallel to face
									continue;
								}
								double lambda = n_pminx / n_xpreviousminx;
								if ((lambda > 1 + 1e-20) or (lambda < -1e-20)) {
									// Lagrangian path does not cut face
									continue;
								}

								dealii::Tensor<1, dim> h = x_previous - x_i;
								h *= lambda;
								dealii::Point<dim> x_b = x_i + h;
								if (not cell->point_inside(x_b)) {
									// another face is cut before this one
									continue;
								}
								found_across_bound = true;
								// Periodic boundaries

								// Other boundaries

							}
						}
					} /* while not bound_found */
//					if (not found_across_bound) {
//						LOG(WARNING) << "The population f" << alpha + 1
//								<< "on support point could not trace its Lagrangian path from the point "
//								<< x_i
//								<< "backwards. Support point not found. As a last resort, all locally "
//										"owned cells are searched recursively. Note that in this way, boundary conditions "
//										"could be omitted. If this is unsuccessful, NATriuM will terminate."
//								<< endl;
//						cout << "The population f" << alpha + 1
//								<< "on support point could not trace its Lagrangian path from the point "
//								<< x_i
//								<< "backwards. Support point not found. As a last resort, all locally "
//										"owned cells are searched recursively. Note that in this way, boundary conditions "
//										"could be omitted. If this is unsuccessful, NATriuM will terminate."
//								<< endl;
//						typename dealii::DoFHandler<dim>::active_cell_iterator found_in_cell =
//								recursivelySearchInNeighborhood(x_previous,
//										cell);
//						if (found_in_cell != m_doFHandler->end()) {
//							// coarse check if boundaries conditions were violated
//							// TODO
//							// add to cell_map
//							DoFInfo info_i_b(local_dof_indices.at(i), alpha,
//									alpha, x_previous);
//							typename std::map<
//									typename DoFHandler<dim>::active_cell_iterator,
//									XList>::iterator it = cell_map.find(
//									found_in_cell);
//							if (it == cell_map.end()) {
//								XList new_xlist;
//								new_xlist.push_back(info_i_b);
//								cell_map.insert(
//										std::pair<
//												typename DoFHandler<dim>::active_cell_iterator,
//												XList>(found_in_cell,
//												new_xlist));
//							} else {
//								it->second.push_back(info_i_b);
//							}
//							break;
//						} else { /* if not found in cell */
//							cout << x_previous << "not found anywhere in the domain." << endl;
//							natrium_errorexit("Semi-Lagrangian assembly failed. A point could not be "
//									"found in the domain. Try to decrease the time step size.");
//						} /* if found in cell // else  */
//					}

				}
			}
		}

	}
	m_systemVector.compress(dealii::VectorOperation::add);
	m_systemMatrix.compress(dealii::VectorOperation::add);

//#endif
} /* fillSparseObject */

//template<>
//int SemiLagrangian<2>::faceCrossedFirst(
//		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
//		const dealii::Point<2>& p_inside, const dealii::Point<2>& p_outside,
//		dealii::Point<2>& p_boundary, double* lambda) {
//	// transform to unit cell
//	dealii::DoFHandler<2>::cell_iterator ci(*cell);
//	dealii::Point<2> pi_unit = m_mapping.transform_real_to_unit_cell(ci,
//			p_inside);
//	dealii::Point<2> po_unit = m_mapping.transform_real_to_unit_cell(ci,
//			p_outside);
//
//	/*       3
//	 *    2-->--3
//	 *    |     |
//	 *   0^     ^1
//	 *    |     |
//	 *    0-->--1
//	 *        2
//	 */
//	int face_id = -1;
//	if (po_unit[0] < 0) {
//		// if face 0 is crossed
//		if (po_unit[1] < 0) {
//			// if face 2 is crossed
//			if (pi_unit[1] / pi_unit[0] > po_unit[1] / po_unit[0])
//				face_id = 2;
//			else
//				face_id = 0;
//		} else if (po_unit[1] > 1) {
//			// if face 3 is crossed
//			if ((1 - po_unit[1]) / po_unit[0] > (1 - pi_unit[1]) / pi_unit[0])
//				face_id = 3;
//			else
//				face_id = 0;
//		} else {
//			face_id = 0;
//		}
//	} else if (po_unit[0] > 1) {
//		// check if face 1 is crossed
//		if (po_unit[1] < 0) {
//			// if face 2 is crossed
//			if (po_unit[1] / (1 - po_unit[0]) > pi_unit[1] / (1 - pi_unit[0]))
//				face_id = 2;
//			else
//				face_id = 1;
//		} else if (po_unit[1] > 1) {
//			// if face 3 is crossed
//			if ((1 - po_unit[1]) / (1 - po_unit[0])
//					> (1 - pi_unit[1]) / (1 - pi_unit[0]))
//				face_id = 3;
//			else
//				face_id = 1;
//
//		} else {
//			// only face 1 is crossed
//			face_id = 1;
//		}
//	} else {
//		if (po_unit[1] < 0) {
//			// only face 2 is crossed
//			face_id = 2;
//		} else if (po_unit[1] > 1) {
//			// only face 3 is crossed
//			face_id = 3;
//		} else {
//			// no face crossed;
//			return -1;
//		}
//	}
//
//	// calculate boundary point
//	dealii::Point<2> pb_unit;
//	size_t coord_b; // boundary coordinate of unit cell (0 or 1)
//	size_t index_b; // index of boundary coordinate ( 0 for x or 1 for 1 )
//	switch (face_id) {
//	case 0: {
//		index_b = 0;
//		coord_b = 0;
//		break;
//	}
//	case 1: {
//		index_b = 0;
//		coord_b = 1;
//		break;
//	}
//	case 2: {
//		index_b = 1;
//		coord_b = 0;
//		break;
//	}
//	case 3: {
//		index_b = 1;
//		coord_b = 1;
//		break;
//	}
//	} /* switch face_id */
//	// calculate lambda and coordinates of boundary point
//	pb_unit[index_b] = coord_b;
//	*lambda = (coord_b - po_unit[coord_b])
//			/ (pi_unit[1 - coord_b] - po_unit[1 - coord_b]);
//	pb_unit[1 - coord_b] = po_unit[1 - coord_b]
//			+ *lambda * (pi_unit[1 - coord_b] - po_unit[1 - coord_b]);
//	p_boundary = m_mapping.transform_unit_to_real_cell(ci, pb_unit);
//
//	assert(cell->point_inside(p_boundary));
//	return face_id;
//} /* faceCrossedFirst */

template<size_t dim>
int SemiLagrangian<dim>::faceCrossedFirst(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const dealii::Point<dim>& p_inside, const dealii::Point<dim>& p_outside,
		dealii::Point<dim>& p_boundary, double* lambda) {
	// transform to unit cell
	typename dealii::DoFHandler<dim>::cell_iterator ci(*cell);
	dealii::Point<dim> pi_unit = m_mapping.transform_real_to_unit_cell(ci,
			p_inside);
	dealii::Point<dim> po_unit = m_mapping.transform_real_to_unit_cell(ci,
			p_outside);

	/*       3
	 *    2-->--3
	 *    |     |
	 *   0^     ^1
	 *    |     |
	 *    0-->--1
	 *        2
	 */

	/*       *-------*        *-------*
	 *      /|       |       /       /|
	 *     / |   3   |      /   5   / |
	 *    /  |       |     /       /  |
	 *   *   |       |    *-------*   |
	 *   | 0 *-------*    |       | 1 *
	 *   |  /       /     |       |  /
	 *   | /   4   /      |   2   | /
	 *   |/       /       |       |/
	 *   *-------*        *-------*
	 *
	 */
	int face_id = -1;
	*lambda = 100;

	// for efficiency reasons: omit the lambda-stuff, if possible (which is in most calls)
	if (po_unit[0] < 0) {
		// if face 0 is crossed
		*lambda = (0 - pi_unit[0]) / (po_unit[0] - pi_unit[0]);
		assert(*lambda > 0);
		assert(*lambda < 1);
		face_id = 0;
	} else if (po_unit[0] > 1) {
		// if face 1 is crossed
		*lambda = (1 - pi_unit[0]) / (po_unit[0] - pi_unit[0]);
		assert(*lambda > 0);
		assert(*lambda < 1);
		face_id = 1;
	}
	if (po_unit[1] < 0) {
		// if face 2 is crossed
		double lambda_y = (0 - pi_unit[1]) / (po_unit[1] - pi_unit[1]);
		assert(lambda_y > 0);
		assert(lambda_y < 1);
		if (lambda_y < *lambda) {
			*lambda = lambda_y;
			face_id = 2;
		}
	} else if (po_unit[1] > 1) {
		// if face 3 is crossed
		double lambda_y = (1 - pi_unit[1]) / (po_unit[1] - pi_unit[1]);
		assert(lambda_y > 0);
		assert(lambda_y < 1);
		if (lambda_y < *lambda) {
			*lambda = lambda_y;
			face_id = 3;
		}
	}
	if (dim == 3) {
		if (po_unit[2] < 0) {
			// if face 4 is crossed
			double lambda_z = (0 - pi_unit[2]) / (po_unit[2] - pi_unit[2]);
			assert(lambda_z > 0);
			assert(lambda_z < 1);
			if (lambda_z < *lambda) {
				*lambda = lambda_z;
				face_id = 4;
			}
		} else if (po_unit[2] > 1) {
			// if face 5 is crossed
			double lambda_z = (1 - pi_unit[2]) / (po_unit[2] - pi_unit[2]);
			assert(lambda_z > 0);
			assert(lambda_z < 1);
			if (lambda_z < *lambda) {
				*lambda = lambda_z;
				face_id = 5;
			}
		}
	}
	if (face_id == -1) {
		// po is inside
		return -1;
	}

	// calculate boundary point
	dealii::Tensor<1, dim> increment = po_unit - pi_unit;
	increment *= (*lambda);
	p_boundary = m_mapping.transform_unit_to_real_cell(ci, pi_unit + increment);

	assert(cell->point_inside(p_boundary));

	// TODO assert that p_boundary is at cell boundary

	return face_id;

} /* faceCrossedFirst */

template class SemiLagrangian<2> ;
template class SemiLagrangian<3> ;

} /* namespace natrium */
