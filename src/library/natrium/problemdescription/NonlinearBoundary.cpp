/*
 * NonlinearBoundary.cpp
 *
 *  Created on: 26.03.2014
 *      Author: kraemer
 */

#include "NonlinearBoundary.h"

#include "deal.II/lac/compressed_sparsity_pattern.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>

namespace natrium {

template<size_t dim> NonlinearBoundary<dim>::NonlinearBoundary(
		size_t boundaryIndicator, const dealii::UpdateFlags update_flags) :
		m_boundaryIndicator(boundaryIndicator), m_updateFlags(update_flags), m_advectionOperator(
		NULL), m_rho(NULL), m_u(NULL), m_f(NULL), m_boundaryVector(NULL) {
}

template<size_t dim>
void NonlinearBoundary<dim>::initialize(
		boost::shared_ptr<AdvectionOperator<dim> > advection_operator,
		boost::shared_ptr<Stencil> stencil,
		distributed_vector const * rho, vector<distributed_vector> const* u,
		DistributionFunctions const * f, distributed_block_vector* boundary_vector) {
	m_advectionOperator = advection_operator;
	m_stencil = stencil;
	m_rho = rho;
	m_u = u;
	m_f = f;
	m_boundaryVector = boundary_vector;
}

template class NonlinearBoundary<2> ;
template class NonlinearBoundary<3> ;

} /* namespace natrium */
