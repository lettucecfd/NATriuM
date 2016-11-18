/*
 * SemiLagrangianBoundaryHandler.cpp
 *
 *  Created on: 15.11.2016
 *      Author: akraem3m
 */

#include "SemiLagrangianBoundaryHandler.h"

namespace natrium {

template<size_t dim>
void SemiLagrangianBoundaryHandler<dim>::addHit(
		const LagrangianPathTracker<dim>& tracker, size_t boundary_id,
		const AdvectionOperator<dim>& sl) {
	assert (m_timeStep > 0);
	assert(m_boundaries.hasID(boundary_id));
	assert(m_boundaries.isSL(boundary_id));
	// check if cell is already there
	typename HitList<dim>::iterator it;
	it = m_hitList.find(tracker.currentCell);
	if (m_hitList.end() == it) {
		HitListAtCell<dim> hl_c;
		typename dealii::DoFHandler<dim>::active_cell_iterator c(
				tracker.currentCell);
		it = m_hitList.insert(std::make_pair(c, hl_c)).first;
	}

	// check if current point is support point
	int q = supportPointNr<dim>(tracker.currentCell, tracker.currentPoint,
			*sl.getQuadrature(), sl.getMapping());
	if (q >= 0) {
		// could be generalized to 'if currentPoint is support point'
		it->second.addHit(
				BoundaryHitAtSupportPoint<dim>(tracker, m_stencil, boundary_id,
						q));
	} else {
		it->second.addHit(BoundaryHit<dim>(tracker, m_stencil, boundary_id));
	}

}

template<size_t dim>
void SemiLagrangianBoundaryHandler<dim>::apply(DistributionFunctions& f_new,
		const DistributionFunctions& f_old, const dealii::DoFHandler<dim>& dof,
		boost::shared_ptr<dealii::FEValues<dim> > fe_values) {

	std::vector<dealii::Point<dim> > off_sup_points;
	std::vector< dealii::types::global_dof_index > local_dofs;
	local_dofs.resize(fe_values->get_fe().n_dofs_per_cell());
	typename HitList<dim>::iterator it = m_hitList.begin();
	typename HitList<dim>::iterator end = m_hitList.end();
	for (; it != end; ++it) {

		const typename dealii::DoFHandler<dim>::active_cell_iterator & cell =
				it->first;
		const HitListAtCell<dim> & cell_hits = it->second;

		// update both fevalues instances, the one for the support points and the one for the arbitrary points
		fe_values->reinit(cell);
		cell->get_dof_indices(local_dofs);

		off_sup_points.clear();
		for (size_t i = 0; i < cell_hits.hitListArbitrary.size(); i++) {
			off_sup_points.push_back(
					cell_hits.hitListArbitrary.at(i).getCurrentPoint());
		}
		boost::shared_ptr<dealii::FEValues<dim> > fe_values2 =
				reinitArbitraryPoints<dim>(cell, off_sup_points,
						fe_values->get_mapping(),
						fe_values->get_update_flags());


		// apply boundaries for support points
		for (size_t i = 0; i < cell_hits.hitListSupportPoints.size(); i++) {
			m_boundaries.getSLBoundary(
					cell_hits.hitListSupportPoints.at(i).getBoundaryId())->calculateBoundaryValues(
					f_old, f_new, local_dofs, *fe_values,
					cell_hits.hitListSupportPoints.at(i).getSupportQPoint(),
					cell_hits.hitListSupportPoints.at(i).getDestination(),
					cell_hits.hitListSupportPoints.at(i).getDtHit());
		}

		// apply boundaries for arbitrary points
		// quadrature is already in the right order by construction, i.e. q==i
		for (size_t q = 0; q < cell_hits.hitListArbitrary.size(); q++) {
			m_boundaries.getSLBoundary(
					cell_hits.hitListSupportPoints.at(q).getBoundaryId())->calculateBoundaryValues(
					f_old, f_new, local_dofs, *fe_values2, q,
					cell_hits.hitListSupportPoints.at(q).getDestination(),
					cell_hits.hitListSupportPoints.at(q).getDtHit());

		}
	}
}

template class BoundaryHit<2> ;
template class BoundaryHit<3> ;
template class BoundaryHitAtSupportPoint<2> ;
template class BoundaryHitAtSupportPoint<3> ;
template class HitListAtCell<2> ;
template class HitListAtCell<3> ;
template class SemiLagrangianBoundaryHandler<2> ;
template class SemiLagrangianBoundaryHandler<3> ;

} /* namespace natrium */
