/*
 * SemiLagrangianBoundaryHandler.cpp
 *
 *  Created on: 15.11.2016
 *      Author: akraem3m
 */

#include "SemiLagrangianBoundaryHandler.h"
#include "BoundaryFlags.h"

namespace natrium {

template<size_t dim>
void SemiLagrangianBoundaryHandler<dim>::addHit(
		const LagrangianPathTracker<dim>& tracker, size_t boundary_id) {
	//pout << m_timeStep << endl;
	assert(m_timeStep > 0);
	assert(m_boundaries.hasID(boundary_id));
	assert(m_boundaries.getBoundary(boundary_id)->isSLSupported());
	// check if cell is already there
	m_hitList.addHit(BoundaryHit<dim>(tracker, m_stencil, boundary_id, m_timeStep));
	// TODO deal with corner nodes. Each boundary hit could be given prescribed values of multiple boundaries.
}

template<size_t dim>
void SemiLagrangianBoundaryHandler<dim>::apply(DistributionFunctions& f_new,
		const DistributionFunctions& f_old, double t) {

	GlobalBoundaryData g_data(f_old, f_new, m_stencil, 0.1, m_timeStep);
	std::vector<dealii::Point<dim> > local_hit_points;
	typename HitList<dim>::iterator cell_it = m_hitList.begin();
	typename HitList<dim>::iterator cell_end = m_hitList.end();
	for (; cell_it != cell_end; ++cell_it) {
		// i.e., FOR ALL CELLS

		const typename dealii::DoFHandler<dim>::active_cell_iterator & cell =
				cell_it->first;
		const HitListAtCell<dim> & cell_hits = cell_it->second;

		// get update flags for all boundaries that are present at this cell
		BoundaryFlags flags = only_distributions;
		for (size_t f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++) {
			if (not cell->at_boundary(f)) {
				continue;
			}
			size_t bi = cell->face(f)->boundary_id();
			if (m_boundaries.getBoundary(bi)->isPeriodic())
				continue;
			if (m_boundaries.getBoundary(bi)->isSLSupported()) {
				BoundaryFlags flags_f = m_boundaries.getBoundary(bi)->getUpdateFlags();
				flags |= flags_f;
			} else {
				throw NATriuMException("SemiLagrangian streaming can only be used with periodic boundaries "
						"and SL boundaries. Your boundary had isSLSupported() == false.");
			}
		}

		// get hit points
		local_hit_points.clear();
		typename HitListAtCell<dim>::const_iterator point_it =
				cell_hits.begin();
		typename HitListAtCell<dim>::const_iterator point_end =
				cell_hits.end();
		for (; point_it != point_end; ++point_it) {
			local_hit_points.push_back(point_it->first);
		}

		// make boundary values instance for the calculations at the boundary
		FEBoundaryValues<dim> fe_b_values(cell,
				local_hit_points,
				*m_advection.getFe(), m_advection.getMapping(),  flags,
				g_data);

		// apply boundaries at boundaries
		size_t q_point = 0;
		for (point_it = cell_hits.begin(); point_it != point_end;
				++point_it) {
			const HitListAtPoint<dim> & point_hits = point_it->second;
			// calculate values at points
			fe_b_values.reinit(q_point);

			// apply boundaries at hits
			for (size_t i = 0; i < point_hits.n_hits(); i++){
				const BoundaryHit<dim>& hit = point_hits.at(i);
				// calculate new distribution functions at hits
				assert (m_boundaries.getBoundary(hit.getBoundaryId())->isSLSupported());
				 m_boundaries.getBoundary(hit.getBoundaryId())->calculateBoundaryValues( fe_b_values,
								q_point, hit.getDestination(),
								hit.getDtHit(), t);
			} /* for all hits */
			q_point ++;
		} /* for all points */
	} /* for all cells */

}

template class SemiLagrangianBoundaryHandler<2> ;
template class SemiLagrangianBoundaryHandler<3> ;

} /* namespace natrium */
