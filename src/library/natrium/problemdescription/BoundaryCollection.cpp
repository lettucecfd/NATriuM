/*
 * BoundaryCollection.cpp
 *
 *  Created on: 29.10.2013
 *      Author: kraemer
 */

#include "BoundaryCollection.h"

namespace natrium {

template<size_t dim>
void BoundaryCollection<dim>::updateNonlinearBoundaryValues() {
	if (!hasNonlinearBoundaries()){
		return;
	}
	BoundaryCollection<dim>::NonlinearIterator it;
	BoundaryCollection<dim>::NonlinearIterator end =
			m_nonlinearBoundaries.end();
	// reset boundary vector to constant (linear) boundary vector
	m_nonlinearBoundaries.begin()->second->resetBoundaryVector();
	for (it = m_nonlinearBoundaries.begin(); it != end; it++) {
		it->second->updateNonlinearBoundaryValues();
	}
}

template<size_t dim>
void BoundaryCollection<dim>::initializeNonlinearBoundaries(
		boost::shared_ptr<AdvectionOperator<dim> > advection_operator,
		boost::shared_ptr<Stencil> stencil, distributed_vector const * rho,
		vector<distributed_vector> const* u, DistributionFunctions const * f,
		distributed_block_vector* boundary_vector) {
	BoundaryCollection<dim>::NonlinearIterator it;
	BoundaryCollection<dim>::NonlinearIterator end =
			m_nonlinearBoundaries.end();
	for (it = m_nonlinearBoundaries.begin(); it != end; it++) {
		it->second->initialize(advection_operator, stencil, rho, u, f,
				boundary_vector);
	}
}

template class BoundaryCollection<2> ;
template class BoundaryCollection<3> ;

} /* namespace natrium */
