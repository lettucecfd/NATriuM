/**
 * @file SemiLagrangian.cpp
 * @short Semi-Lagrangian advection scheme
 * @date 29.04.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

//#include "SLBoundary.h"
#include "SemiLagrangian.h"

// #include <fstream>
#include <queue>
#include <array>
#include <sstream>

#include "deal.II/lac/trilinos_sparsity_pattern.h"
#include "deal.II/dofs/dof_renumbering.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/lac/matrix_iterator.h"
#include "deal.II/lac/sparsity_tools.h"
#include "deal.II/base/utilities.h"
#include "deal.II/fe/fe_q.h"
#include "deal.II/fe/mapping_q1.h"
#include "deal.II/grid/grid_tools.h"

#include "../boundaries/PeriodicBoundary.h"
#include "../stencils/Stencil.h"

#include "../utilities/DealiiExtensions.h"
#include "../utilities/CFDSolverUtilities.h"

using namespace dealii;

namespace natrium {

template<size_t dim>
SemiLagrangian<dim>::SemiLagrangian(ProblemDescription<dim>& problem,
		size_t orderOfFiniteElement, QuadratureName quad_name,
		SupportPointsName points_name, boost::shared_ptr<Stencil> stencil,
		double delta_t) :
		AdvectionOperator<dim>(problem, orderOfFiniteElement, quad_name,
				points_name, stencil, delta_t, false), m_boundaryHandler(*this) {

	m_boundaryHandler.setTimeStep(delta_t);

} /* SEDGMinLee<dim>::SEDGMinLee */

template<size_t dim>
SemiLagrangian<dim>::SemiLagrangian(ProblemDescription<dim>& problem,
		size_t orderOfFiniteElement, boost::shared_ptr<Stencil> stencil,
		double delta_t) :
		SemiLagrangian(problem, orderOfFiniteElement, QGAUSS_LOBATTO,
				GAUSS_LOBATTO_POINTS, stencil, delta_t) {

	m_boundaryHandler.setTimeStep(delta_t);

}

template<size_t dim>
void SemiLagrangian<dim>::setupDoFs() {

	LOG(DETAILED) << "Setup DoFs..." << endl;
// distribute degrees of freedom over mesh
	Base::distributeDoFs();

	updateSparsityPattern();

	LOG(DETAILED) << "... done (setup DoFs)." << endl;

}

template<size_t dim>
void SemiLagrangian<dim>::reassemble() {
// TODO: if Mesh changed: reinit dof-handler and sparsity pattern in some way

	LOG(DETAILED) << "Assemble..." << endl;
// make sure that sparsity structure is not empty
	assert(Base::m_systemMatrix.n() != 0);
	assert(Base::m_systemMatrix.m() != 0);

	fillSparseObject(false); // "false" says that you want to fill the system matrix, not the sparsity pattern

	Base::m_systemMatrix.compress(dealii::VectorOperation::add);

	LOG(DETAILED) << "... done (assemble)." << endl;

} /* reassemble */

template<size_t dim>
void SemiLagrangian<dim>::updateSparsityPattern() {

	LOG(DETAILED) << "Update sparsity pattern..." << endl;
// TODO only update sparsity pattern for changed cells
///////////////////////////////////////////////////////
// Setup sparsity pattern (completely manually):
///////////////////////////////////////////////////////
// allocate sizes
	size_t n_blocks = Base::m_stencil->getQ() - 1;

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
			m_sparsityPattern[i][j].reinit(Base::getLocallyOwnedDofs(),
					Base::getLocallyOwnedDofs(), Base::getLocallyRelevantDofs(),
					MPI_COMM_WORLD);

// fill sparsity pattern
	if (Base::m_deltaT != 0) {
		fillSparseObject(true); // "true" means that the object to fill is the sparsity pattern
	}

// initialize matrix
	Base::m_systemMatrix.reinit(n_blocks, n_blocks);
	for (size_t i = 0; i < n_blocks; i++) {
		for (size_t j = 0; j < n_blocks; j++) {
			m_sparsityPattern[i][j].compress();
			Base::m_systemMatrix.block(i, j).reinit(m_sparsityPattern[i][j]);
		}
	}
	Base::m_systemMatrix.collect_sizes();

	LOG(DETAILED) << "... done (update sparsity pattern)." << endl;
}
/* updateSparsityPattern */

template<size_t dim>
void SemiLagrangian<dim>::fillSparseObject(bool sparsity_pattern) {

//TimerOutput::Scope timer_section(Timing::getTimer(),
//			"Assembly: fill sparse object");

	dealii::IndexSet tracked_indices(Base::m_doFHandler->n_dofs());
	if (not sparsity_pattern) {
		// make sure that sparsity structure is not empty
		assert(Base::m_systemMatrix.n() != 0);
		assert(Base::m_systemMatrix.m() != 0);

	}
	assert(m_boundaryHandler.n_hits() == 0);

/////////////////////////////////
// Initialize Finite Element ////

// Initialize
	std::map<typename DoFHandler<dim>::active_cell_iterator,
			DeparturePointList<dim> > found_in_cell;
	const size_t dofs_per_cell = Base::m_fe->dofs_per_cell;
	const std::vector<Point<dim> > & unit_support_points =
			Base::m_fe->get_unit_support_points();
	std::vector<std::vector<double> > local_entries;
	std::vector<Point<dim> > local_lagrange_points;
	std::vector<Tensor<1, dim> > minus_dtealpha;
	for (size_t i = 0; i < Base::m_stencil->getQ(); i++) {
		minus_dtealpha.push_back(
				vectorToTensor<dim>(Base::m_stencil->getDirection(i))
						* (-Base::m_deltaT));
	}
	std::vector<dealii::types::global_dof_index> local_dof_indices(
			dofs_per_cell);
	std::vector<dealii::types::global_dof_index> neighbor_dof_indices(
			dofs_per_cell);
	size_t max_n_shells =
			(size_t) (Base::m_stencil->getMaxParticleVelocityMagnitude()
					* Base::m_deltaT
					/ CFDSolverUtilities::getMinimumVertexDistance<dim>(
							*Base::m_problem.getMesh())) + 1;
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
			Base::m_doFHandler->begin_active(), endc =
			Base::m_doFHandler->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// initialize
			found_in_cell.clear();
			std::queue<LagrangianPathTracker<dim> > not_found;

			// get global degrees of freedom
			cell->get_dof_indices(local_dof_indices);


			// ================================================================================================
			// =================================  Create Lagrangian support points ============================
			// ================================================================================================
			// for all support points in cell
			for (size_t i = 0; i < dofs_per_cell; i++) {
				// force each support point to be handled only once
				if (tracked_indices.is_element(local_dof_indices.at(i))) {
					continue;
				}
				if (not Base::m_doFHandler->locally_owned_dofs().is_element(
						local_dof_indices.at(i))) {
					continue;
				}
				tracked_indices.add_index(local_dof_indices.at(i));
				//TimerOutput::Scope timer_section(Timing::getTimer(),
				//			"Assembly: create points");
				// get a point x
				dealii::Point<dim> x_i =
						StaticMappingQ1<dim,dim>::mapping.transform_unit_to_real_cell(cell,
								unit_support_points.at(i));
				// for all directions
				for (size_t alpha = 1; alpha < Base::m_stencil->getQ();
						alpha++) {
					// calculate x^(t-delta_t)
					dealii::Point<dim> x_departure = x_i
							+ minus_dtealpha.at(alpha);
					LagrangianPathTracker<dim> info_i(local_dof_indices.at(i),
							alpha, alpha, x_departure, x_i, cell);
					not_found.push(info_i);
				}
			}

			// -- Recursively follow Lagrangian paths ---
			// while there are points that have not been found, yet
			while (not not_found.empty()) {
				//TimerOutput::Scope timer_section(Timing::getTimer(),
				//			"Assembly: follow paths");
				LagrangianPathTracker<dim>& el = not_found.front();
				el.lifeTimeCounter++;
				// check if el is running into Nirvana (which can happen, e.g., at domain corners)
				if (el.lifeTimeCounter > 50)
				{
					// write '1' on diagonal (i.e., do not change this distribution function during streaming step)
					if (sparsity_pattern) {
						// add diagonal entry to sparsity pattern
						m_sparsityPattern[el.destination.direction - 1][el.destination.direction - 1].add(el.destination.index,
								el.destination.index);
					} else {
						// insert '1'
						Base::m_systemMatrix.block(
								el.destination.direction - 1,
								el.destination.direction - 1).add(
								el.destination.index, el.destination.index,
								1.0);
					}
					not_found.pop();
					continue;
				}

				// ================================================================================================
				// ================================= Find out which (and if) face is crossed ======================
				// ================================================================================================
				// find out which face is crossed first
				double lambda = 100;
				size_t child_id = 100;
				dealii::Point<dim> p_boundary;
				int face_id = -2;
				try{
					face_id = faceCrossedFirst(el.currentCell, el.currentPoint,
							el.departurePoint, p_boundary, &lambda, &child_id);
				} catch (FaceCrossedFirstFailed& e){
					try {
						el.currentCell = dealii::GridTools::find_active_cell_around_point(*(Base::m_problem.getMesh()), 	el.departurePoint);
						el.currentPoint = el.departurePoint;
						continue;
					} catch( dealii::GridTools::ExcPointNotFound<dim>& f) {
						std::stringstream msg;
						msg << "NATriuM had to rely on deal.II's function find_active_cell_around point "
								" because NATriuM's function face_crossed_first did not work due to a non-invertible "
								" mapping function. The deal.II function found out that the point does not "
								" lie inside the grid. This situation is not supported by the current version of the Semi-Lagrangian code. "
								" TODO: implement a supplement for face_crossed_first that does not rely on inverting the mapping function."
								<< endl;
						LOG(ERROR) << msg;
						throw SemiLagrangianException(msg.str());
					}
				}
				if (face_id == -1) {
					// point found in this cell: add to cell_map
					typename std::map<
							typename DoFHandler<dim>::active_cell_iterator,
							DeparturePointList<dim> >::iterator it =
							found_in_cell.find(el.currentCell);
					if (it == found_in_cell.end()) {
						DeparturePointList<dim> new_xlist;
						new_xlist.push_back(el);
						found_in_cell.insert(
								std::pair<
										typename DoFHandler<dim>::active_cell_iterator,
										DeparturePointList<dim> >(
										el.currentCell, new_xlist));
					} else {
						it->second.push_back(el);
					}
					not_found.pop();

				// ================================================================================================
				// ================================= Boundary faces ===============================================
				// ================================================================================================
				} else if (el.currentCell->at_boundary(face_id)) {
					size_t bi = el.currentCell->face(face_id)->boundary_id();
					if (Base::getBoundaries()->isPeriodic(bi)) {
						// Apply periodic boundaries
						const boost::shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
								Base::getBoundaries()->getPeriodicBoundary(bi);
						assert(
								periodicBoundary->isFaceInBoundary(
										el.currentCell, face_id));
						el.currentPoint =
								periodicBoundary->coordinatesAcrossPeriodicBoundary(
										p_boundary, el.currentCell);
						el.departurePoint =
								periodicBoundary->coordinatesAcrossPeriodicBoundary(
										el.departurePoint, el.currentCell);
						const typename dealii::DoFHandler<dim>::active_cell_iterator h =
								el.currentCell;
						periodicBoundary->getOppositeCellAtPeriodicBoundary(h,
								el.currentCell);

					} else /* if is not periodic */{
						if ((Base::getBoundaries()->getBoundary(bi)->isLinearFluxBoundary())) {

							for (size_t i = 0; i < dim; i++) {
															pout << "curr" << i << " " << el.currentPoint[i]<< endl;
															pout << "p_b" << i << " " << p_boundary[i]<< endl;
															pout << "Dep "<< i << " = " << el.departurePoint[i] << endl;
														}

							pout << "Direction "<< el.currentDirection << endl;

							el.currentDirection =
									Base::m_stencil->getIndexOfOppositeDirection(
											el.currentDirection);
							el.currentPoint = p_boundary;

							double vel_direction = 0;
							for (size_t i = 0; i < dim; i++) {
								vel_direction += Base::m_stencil->getDirection(
										el.currentDirection)[i]
										* Base::m_stencil->getDirection(
												el.currentDirection)[i];
							}
							vel_direction = sqrt(vel_direction);
							
							double distance = el.departurePoint.distance(
									el.currentPoint);
							for (size_t i = 0; i < dim; i++) {
								el.departurePoint[i] = el.currentPoint[i]
										- Base::m_stencil->getDirection(
												el.currentDirection)[i]
												* distance / vel_direction;
							}

							for (size_t i = 0; i < dim; i++) {
																						pout << "Depa"<< i << " = " << el.departurePoint[i] << endl;
																					}
							pout << endl;

						//	m_boundaryHandler.addHit(el, bi);
							// else
							if (not sparsity_pattern) {
								el.currentPoint = p_boundary;
							//	m_boundaryHandler.addHit(el, bi);
							}   // not_found.pop();

						} //endif isLinearFluxBoundary

					} /* endif isPeriodic */
				// ================================================================================================
				// ================================= Interior faces ===============================================
				// ================================================================================================
				} else {
					// Interior faces
					el.currentPoint = p_boundary;
					if (el.currentCell->neighbor(face_id)->has_children()) {
						el.currentCell =
								el.currentCell->neighbor(face_id)->child(
										child_id);
					} else {
						el.currentCell = el.currentCell->neighbor(face_id);
					}
					if (not (el.currentCell->is_locally_owned()
							or el.currentCell->is_ghost())) {
						// information at source point is on a different processor
						std::stringstream s;
						s
								<< "Time step too large in semi-Lagrangian streaming. Lagrangian source point "
								<< el.departurePoint
								<< " is owned by another MPI process, but required by MPI process "
								<< dealii::Utilities::MPI::this_mpi_process(
								MPI_COMM_WORLD)
								<< ". Please decrease the CFL number so that the Lagrangian source"
								<< " points are at least inside the layer of ghost cells."
								<< endl;
						natrium_errorexit(s.str().c_str());
					} /* end if cell not locally owned or ghost */
				} /* end if at boundary / else */

			} /* end while not found */

			// ================================================================================================
			// ================================= Assemble matrix/sparsity pattern =============================
			// ================================================================================================
			typename std::map<typename DoFHandler<dim>::active_cell_iterator,
					DeparturePointList<dim> >::iterator list =
					found_in_cell.begin();
			typename std::map<typename DoFHandler<dim>::active_cell_iterator,
					DeparturePointList<dim> >::iterator end_list =
					found_in_cell.end();
			// (i.e. for all cells that contain support points )
			for (; list != end_list; list++) {
				//TimerOutput::Scope timer_section(Timing::getTimer(),
				//			"Assembly: add entries");
				const typename DoFHandler<dim>::active_cell_iterator & it =
						list->first;

				// get dofs
				it->get_dof_indices(local_dof_indices);
				const DeparturePointList<dim> & l = list->second;

				//if (not sparsity_pattern) {
				// calculate shape values
				local_entries.clear();
				local_lagrange_points.clear();
				for (size_t i = 0; i < l.size(); i++) {
					local_lagrange_points.push_back(l[i].departurePoint);
					std::vector<double> h(dofs_per_cell);
					local_entries.push_back(h);
				}
				shapeFunctionValue<dim>(it, local_lagrange_points,
						local_entries, *Base::m_mapping);
				//}
				// (i.e. for all Lagrangian points in cell)
				for (size_t i = 0; i < l.size(); i++) {
					assert(it == l.at(i).currentCell);
					for (size_t j = 0; j < local_entries.at(i).size(); j++) {
						if (fabs(local_entries.at(i).at(j)) < 1e-10) {
							continue;
						}
						if (sparsity_pattern) {
							// add entry to sparsity pattern
							m_sparsityPattern[l[i].destination.direction - 1][l[i].currentDirection
									- 1].add(l[i].destination.index,
									local_dof_indices.at(j));
						} else {
							// calculate matrix entries
							//TimerOutput::Scope timer_section(Timing::getTimer(),
							//		"Assembly: block and add");
							Base::m_systemMatrix.block(
									l[i].destination.direction - 1,
									l[i].currentDirection - 1).add(l[i].destination.index,
									local_dof_indices.at(j),
									local_entries.at(i).at(j));	//fe_cell_values.shape_value;
						}
					}
				} /* end for all support points in cell*/
			} /* end for all xlists (i.e. for all cells that contain support points )*/
		} /* end if cell is locally owned */
	} /* end for all cells */

	/*is done later
	 * if (sparsity_pattern) {
	 const size_t Q = Base::m_stencil->getQ();
	 for (size_t i = 0; i < Q - 1; i++) {
	 for (size_t j = 0; j < Q - 1; j++) {
	 m_sparsityPattern.at(i).at(j).compress();
	 }
	 }
	 } else {
	 m_systemVector.compress(dealii::VectorOperation::add);
	 Base::Base::m_systemMatrix.compress(dealii::VectorOperation::add);
	 }*/

//#endif
} /* fillSparseObject */

template<size_t dim>
int SemiLagrangian<dim>::faceCrossedFirst(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const dealii::Point<dim>& p_inside, const dealii::Point<dim>& p_outside,
		dealii::Point<dim>& p_boundary, double* lambda, size_t* child_id) {

//TimerOutput::Scope timer_section(Timing::getTimer(),

// transform to unit cell
	typename dealii::DoFHandler<dim>::cell_iterator ci(*cell);
	dealii::Point<dim> pi_unit = StaticMappingQ1<dim,dim>::mapping.transform_real_to_unit_cell(
			ci, p_inside);
	dealii::Point<dim> po_unit;
	try {
		// the bilinear mapping function might not be invertible outside the cell
		po_unit =
				StaticMappingQ1<dim, dim>::mapping.transform_real_to_unit_cell(
						ci, p_outside);
	} catch (dealii::Mapping<2, 2>::ExcTransformationFailed& e){
		LOG(WARNING)
				<< " using cell_around_active_point is a workaround! "
						"Close to boundaries, the departure point (p_outside) could either reside outside the cell, "
						"or intermediate faces could be missed."
						"Both situations are not handled by this workaround."
				<< endl;
		throw FaceCrossedFirstFailed("Could not calculate inverse mapping.");
	}

// eliminate round-off-errors
	for (size_t i = 0; i < dim; i++) {
		if (fabs(pi_unit[i]) < 1e-10) {
			pi_unit[i] = 0;
		}
		if (fabs(pi_unit[i] - 1) < 1e-10) {
			pi_unit[i] = 1;
		}
		if (fabs(po_unit[i]) < 1e-10) {
			po_unit[i] = 0;
		}
		if (fabs(po_unit[i] - 1) < 1e-10) {
			po_unit[i] = 1;
		}
	}
// make sure the inner point is inside the cell
	for (size_t i = 0; i < dim; i++) {
		if (pi_unit[i] > 1) {
			cout
					<< "The current point does not seem to lie inside the current cell."
							"This is a rather curious error that I encountered once or twice "
							"when working on grids that were not refined at all. I did not spend "
							"too much time to try to make it work. I might include some "
							"assertions in SemiLagrangian.cpp to check where this comes from --"
					<< endl;
			cout << "but not now :-/" << endl;
			cout << "piunit: " << pi_unit[i] << endl;
		}
		assert(pi_unit[i] >= 0);
		assert(pi_unit[i] <= 1);

	}

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
		assert(*lambda >= 0);
		assert(*lambda <= 1);
		face_id = 0;
	} else if (po_unit[0] > 1) {
		// if face 1 is crossed
		*lambda = (1 - pi_unit[0]) / (po_unit[0] - pi_unit[0]);
		assert(*lambda >= 0);
		assert(*lambda <= 1);
		face_id = 1;
	}
	if (po_unit[1] < 0) {
		// if face 2 is crossed
		double lambda_y = (0 - pi_unit[1]) / (po_unit[1] - pi_unit[1]);
		assert(lambda_y >= 0);
		assert(lambda_y <= 1);
		if (lambda_y < *lambda) {
			*lambda = lambda_y;
			face_id = 2;
		}
	} else if (po_unit[1] > 1) {
		// if face 3 is crossed
		double lambda_y = (1 - pi_unit[1]) / (po_unit[1] - pi_unit[1]);
		assert(lambda_y >= 0);
		assert(lambda_y <= 1);
		if (lambda_y < *lambda) {
			*lambda = lambda_y;
			face_id = 3;
		}
	}
	if (dim == 3) {
		if (po_unit[2] < 0) {
			// if face 4 is crossed
			double lambda_z = (0 - pi_unit[2]) / (po_unit[2] - pi_unit[2]);
			assert(lambda_z >= 0);
			assert(lambda_z <= 1);
			if (lambda_z < *lambda) {
				*lambda = lambda_z;
				face_id = 4;
			}
		} else if (po_unit[2] > 1) {
			// if face 5 is crossed
			double lambda_z = (1 - pi_unit[2]) / (po_unit[2] - pi_unit[2]);
			assert(lambda_z >= 0);
			assert(lambda_z <= 1);
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

// compensate for round-off errors
	dealii::Point<dim> h = pi_unit + increment;
	for (size_t i = 0; i < dim; i++) {
		if (fabs(h[i]) < 1e-10) {
			h[i] = 0;
		}
		if (fabs(h[i] - 1) < 1e-10) {
			h[i] = 1;
		}
	}
// assert that point is at boundary
	if (2 == dim) {
		assert(h[0] * h[1] * (1 - h[0]) * (1 - h[1]) == 0);
	} else { // 3 == dim
		assert(h[0] * h[1] * h[2] * (1 - h[0]) * (1 - h[1]) * (1 - h[2]) == 0);
	}

// map to real cell
	p_boundary = StaticMappingQ1<dim,dim>::mapping.transform_unit_to_real_cell(ci, h);

// cout << "boundary point: " << p_boundary << endl;

//assert(cell->point_inside(p_boundary));
	*child_id = dealii::GeometryInfo<dim>::child_cell_from_point(h);

	return face_id;

} /* faceCrossedFirst */

template<size_t dim>
typename dealii::DoFHandler<dim>::cell_iterator SemiLagrangian<dim>::getNeighbor(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t i) {
// cell at periodic boundary
	if (cell->face(i)->at_boundary()) {
		size_t boundaryIndicator = cell->face(i)->boundary_id();
		if (Base::getBoundaries()->isPeriodic(boundaryIndicator)) {
			const boost::shared_ptr<PeriodicBoundary<dim> >& periodicBoundary =
					Base::getBoundaries()->getPeriodicBoundary(
							boundaryIndicator);
			assert(periodicBoundary->isFaceInBoundary(cell, i));
			typename dealii::DoFHandler<dim>::cell_iterator neighborCell;
			periodicBoundary->getOppositeCellAtPeriodicBoundary(cell,
					neighborCell);
			return neighborCell;
		} else {
			return cell->get_dof_handler().end();
		}
	}
// cell not at periodic boundary
	return cell->neighbor(i);
} /* getNeighbor */

/**
 * @short fill the neighborhood list
 * @param cell cell
 * @param neighborhood the neighborhood object
 * @note the neighborhood incorporates the current cell, all its neighbors,
 * 		 and their respective neighbors; each cell has only pointer to it in the neighborhood.
 */
template<size_t dim>
void SemiLagrangian<dim>::getNeighborhood(
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		Neighborhood<dim>& neighborhood, size_t n_shells) {
	const size_t faces_per_cell = 2 * dim;
	neighborhood.clear();
	std::vector < std::array<size_t, faces_per_cell> > visited_faces;

// add self
	neighborhood.push_back(cell);
	visited_faces.push_back(std::array<size_t, faces_per_cell>());
// as the neighborhood gets extended every time a new neighbor is added
// this loop runs over all cells (therefor we have to have a break)
	for (size_t c = 0; c < neighborhood.size(); c++) {
		for (size_t i = 0; i < faces_per_cell; i++) {
			// TODO check the shells over cell orientation rather than face id
			// The face-id test is risky for arbitrary meshes
			if (visited_faces.at(c).at(i) >= n_shells) {
				continue;
			}
			typename dealii::DoFHandler<dim>::cell_iterator n = getNeighbor(
					neighborhood.at(c), i);
			if (n == cell->get_dof_handler().end()) {
				continue;
			} else if (n->is_artificial()) {
				continue;
			} else if (isCellInNeighborhood(*n, neighborhood)) {
				continue;
			} else if (n->active() or n->is_ghost()) {
				std::array < size_t, faces_per_cell > v(visited_faces.at(c));
				v[i]++;
				visited_faces.push_back(v);
				neighborhood.push_back(
						typename dealii::DoFHandler<dim>::active_cell_iterator(
								n));
				continue;
			} else if (n->has_children()) {
				for (size_t j = 0; j < n->n_children(); j++) {
					assert(n->child(j)->active());
					if (isCellInNeighborhood(*(n->child(j)), neighborhood)) {
						continue;
					}
					std::array < size_t, faces_per_cell > v(visited_faces.at(c));
					v[i]++;
					visited_faces.push_back(v);
					neighborhood.push_back(
							typename dealii::DoFHandler<dim>::active_cell_iterator(
									n->child(j)));
					continue;
				}
			}
		}
	}
} /* getNeighborhood */

/**
 * @short recursively search a point in neighborhood, until is found
 * @param p the point you search for
 * @param cell The start cell of the recursive search
 * @return A cell that contains the point p. If the point was not found, the cell pointer will point to DoFHandler.end()
 */
template<size_t dim>
typename dealii::DoFHandler<dim>::active_cell_iterator SemiLagrangian<dim>::recursivelySearchInNeighborhood(
		const dealii::Point<dim>& p,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell) {
	Neighborhood<dim> neighborhood;
// add self
	neighborhood.push_back(cell);
	for (size_t c = 0; c < neighborhood.size(); c++) {
		if (neighborhood.at(c)->point_inside(p)) {
			return neighborhood.at(c);
		}
		for (size_t i = 0; i < dealii::GeometryInfo<dim>::faces_per_cell; i++) {
			typename dealii::DoFHandler<dim>::cell_iterator n = getNeighbor(
					neighborhood.at(c), i);
			if (n == cell->get_dof_handler().end()) {
				continue;
			} else if (n->is_artificial()) {
				continue;
			} else if (isCellInNeighborhood(*n, neighborhood)) {
				continue;
			} else if (n->active() or n->is_ghost()) {
				neighborhood.push_back(
						typename dealii::DoFHandler<dim>::active_cell_iterator(
								n));
				continue;
			} else if (n->has_children()) {
				for (size_t j = 0; j < n->n_children(); j++) {
					assert(n->child(j)->active());
					if (isCellInNeighborhood(*(n->child(j)), neighborhood)) {
						continue;
					}
					neighborhood.push_back(
							typename dealii::DoFHandler<dim>::active_cell_iterator(
									n->child(j)));
					continue;
				}
			}
		}
	}
	return cell->get_dof_handler().end();
} /* recursivelySearchInNeighborhood */

template class SemiLagrangian<2> ;
template class SemiLagrangian<3> ;

} /* namespace natrium */
