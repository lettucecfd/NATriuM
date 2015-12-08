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
		boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
		m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
				boundaryDensity), m_boundaryVelocity(boundaryVelocity) {
}
template DirichletBoundary<2>::DirichletBoundary(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<2> > boundaryDensity,
		boost::shared_ptr<dealii::Function<2> > boundaryVelocity);
template DirichletBoundary<3>::DirichletBoundary(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<3> > boundaryDensity,
		boost::shared_ptr<dealii::Function<3> > boundaryVelocity);

template<size_t dim>
DirichletBoundary<dim>::DirichletBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity) :
		m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
				boost::make_shared<BoundaryDensity<dim> >()), m_boundaryVelocity(
				boost::make_shared<BoundaryVelocity<dim> >(velocity)) {
}
template DirichletBoundary<2>::DirichletBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity);
template DirichletBoundary<3>::DirichletBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity);



} /* namespace natrium */
