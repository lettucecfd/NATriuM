/*
 * CellShockSensor.cpp
 *
 *  Created on: 16.09.2019
 *      Author: natrium
 */

#include "CellShockSensor.h"
#include "deal.II/fe/fe_tools.h"
#include "deal.II/base/polynomial_space.h"
#include "deal.II/base/polynomial.h"

namespace natrium {
template<size_t dim>
CellShockSensor<dim>::CellShockSensor(const dealii::FiniteElement<dim>& fe):m_sourceFE(fe) {
	// TODO Auto-generated constructor stub
}
template<size_t dim>
CellShockSensor<dim>::~CellShockSensor() {
	// TODO Auto-generated destructor stub
}

template<size_t dim>
void CellShockSensor<dim>::applyFilter(
		const dealii::DoFHandler<dim>& dof_handler,
		distributed_vector& dof_vector) {

	const size_t dofs_per_cell = m_sourceFE.dofs_per_cell;
	std::vector<dealii::types::global_dof_index> local_dof_indices(
			dofs_per_cell);
	numeric_vector source_fe_dofs(dofs_per_cell);

	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// get global degrees of freedom
			cell->get_dof_indices(local_dof_indices);

		}
		} /* if is locally owned */

	} /* for cells */


} /* namespace natrium */
