/*
 * DirichletBoundary.cpp
 *
 *  Created on: 26.03.2014
 *      Author: kraemer
 */

#include "deal.II/lac/compressed_sparsity_pattern.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include "DirichletBoundary.h"

namespace natrium {

template class BoundaryDensity<2> ;
template class BoundaryDensity<3> ;
template class BoundaryVelocity<2> ;
template class BoundaryVelocity<3> ;

template<size_t dim> DirichletBoundary<dim>::DirichletBoundary(
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

template<size_t dim>
DirichletBoundary<dim>::DirichletBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity,
		BoundaryTools::DistributionCouplingAtBoundary distribution_coupling,
		BoundaryTools::PointCouplingAtBoundary point_coupling) :
		m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
				boost::make_shared<BoundaryDensity<dim> >()), m_boundaryVelocity(
				boost::make_shared<BoundaryVelocity<dim> >(velocity)),
				m_distributionCoupling(distribution_coupling),
				m_pointCoupling(point_coupling) {
}


template<size_t dim> void DirichletBoundary<dim>::addToSparsityPattern(
		dealii::TrilinosWrappers::SparsityPattern& cSparse,
		const dealii::DoFHandler<dim>& doFHandler) const {
	BoundaryTools::CoupleDoFsAtBoundary<dim>(cSparse,
			doFHandler, m_boundaryIndicator, m_pointCoupling);
}

template class DirichletBoundary<2>;
template class DirichletBoundary<3>;


} /* namespace natrium */
