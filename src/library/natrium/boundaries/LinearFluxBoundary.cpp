/*
 * LinearBoundary.cpp
 *
 *  Created on: 26.03.2014
 *      Author: kraemer
 */

#include "LinearFluxBoundary.h"

#include "deal.II/lac/compressed_sparsity_pattern.h"
#include "deal.II/dofs/dof_handler.h"


namespace natrium {


template<size_t dim> LinearFluxBoundary<dim>::LinearFluxBoundary(
		size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
		boost::shared_ptr<dealii::Function<dim> > boundaryVelocity,
		BoundaryTools::DistributionCouplingAtBoundary distribution_coupling,
		BoundaryTools::PointCouplingAtBoundary point_coupling) :
		m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
				boundaryDensity), m_boundaryVelocity(boundaryVelocity),
				m_distributionCoupling(distribution_coupling),
				m_pointCoupling(point_coupling){
}


template<size_t dim> void LinearFluxBoundary<dim>::addToSparsityPattern(
		dealii::TrilinosWrappers::SparsityPattern& cSparse,
		const dealii::DoFHandler<dim>& doFHandler) const {
	BoundaryTools::CoupleDoFsAtBoundary<dim>(cSparse,
			doFHandler, m_boundaryIndicator, m_pointCoupling);
}

template class LinearFluxBoundary<2>;
template class LinearFluxBoundary<3>;


} /* namespace natrium */
