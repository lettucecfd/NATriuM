/**
 * @file PhysicalProperties.cpp
 * @short
 * @date 06.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "PhysicalProperties.h"

namespace natrium {

template<> PhysicalProperties<2>::PhysicalProperties() {
}
template<> PhysicalProperties<3>::PhysicalProperties() {
}
template<> PhysicalProperties<2>::~PhysicalProperties() {
}
template<> PhysicalProperties<3>::~PhysicalProperties() {
}

template<> double PhysicalProperties<2>::kineticEnergy(
		const vector<distributed_vector>& u, const distributed_vector& rho,
		shared_ptr<AdvectionOperator<2> > advection) {
	const size_t n_dofs = u.at(0).size();
	assert(n_dofs == rho.size());
	assert(n_dofs == u.at(1).size());

	const distributed_vector& ux = u.at(0);
	const distributed_vector& uy = u.at(1);

	// Integrate ux over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values;
	const dealii::DoFHandler<2> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<2> feCellValues(advection->getMapping(),
			*(advection->getFe()), *(advection->getQuadrature()),
			cellUpdateFlags);
	double result = 0.0;
	size_t dofs_per_cell = advection->getFe()->dofs_per_cell;

	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {

			// get global degrees of freedom
			std::vector<dealii::types::global_dof_index> localDoFIndices(
					dofs_per_cell);
			cell->get_dof_indices(localDoFIndices);
			// calculate the fe values for the cell
			feCellValues.reinit(cell);

			size_t local_i;
			for (size_t i = 0; i < dofs_per_cell; i++) {
				local_i = localDoFIndices.at(i);
				result +=
						rho(local_i)
								* (ux(local_i) * ux(local_i)
										+ uy(local_i) * uy(local_i))
								* feCellValues.JxW(
										advection->getCelldofToQIndex().at(i));
			}
		} /* if is locally owned */
	} /* for cells */

	// communicate among MPI processes
	dealii::Utilities::MPI::MinMaxAvg global_res = dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	return 0.5 * global_res.sum;

}
template<> double PhysicalProperties<3>::kineticEnergy(
		const vector<distributed_vector>& u, const distributed_vector& rho,
		shared_ptr<AdvectionOperator<3> > advection) {
	const size_t n_dofs = u.at(0).size();
	assert(n_dofs == rho.size());
	assert(n_dofs == u.at(1).size());
	assert(n_dofs == u.at(2).size());

	const distributed_vector& ux = u.at(0);
	const distributed_vector& uy = u.at(1);
	const distributed_vector& uz = u.at(2);

	// Integrate ux over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values;
	const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<3> feCellValues(advection->getMapping(),
			*(advection->getFe()), *(advection->getQuadrature()),
			cellUpdateFlags);
	double result = 0.0;
	size_t dofs_per_cell = advection->getFe()->dofs_per_cell;

	typename dealii::DoFHandler<3>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// get global degrees of freedom
			std::vector<dealii::types::global_dof_index> localDoFIndices(
					dofs_per_cell);
			cell->get_dof_indices(localDoFIndices);
			// calculate the fe values for the cell
			feCellValues.reinit(cell);

			size_t local_i;
			for (size_t i = 0; i < dofs_per_cell; i++) {
				local_i = localDoFIndices.at(i);
				result += rho(local_i)
						* (ux(local_i) * ux(local_i) + uy(local_i) * uy(local_i)
								+ uz(local_i) * uz(local_i))
						* feCellValues.JxW(
								advection->getCelldofToQIndex().at(i));
			}
		} /* if is locally owned */
	} /* for cells */

	// communicate among MPI processes
	dealii::Utilities::MPI::MinMaxAvg global_res = dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	return 0.5 * global_res.sum;
}

template<size_t dim>
double PhysicalProperties<dim>::maximalPressure(const distributed_vector& rho,
		const double speedOfSound, double & minimalPressure) {
	double maximal = rho.max();
	minimalPressure = rho.min();
	minimalPressure = (minimalPressure - 1.0) * speedOfSound * speedOfSound;
	return (maximal - 1.0) * (speedOfSound * speedOfSound);

}
template
double PhysicalProperties<2>::maximalPressure(const distributed_vector& rho,
		const double speedOfSound, double & minimalPressure);
template
double PhysicalProperties<3>::maximalPressure(const distributed_vector& rho,
		const double speedOfSound, double & minimalPressure);

template<size_t dim>
double PhysicalProperties<dim>::meanVelocityX(const distributed_vector& ux,
		shared_ptr<AdvectionOperator<dim> > advection) {

	// Integrate ux over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values;
	const dealii::DoFHandler<dim> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<dim> feCellValues(advection->getMapping(),
			*(advection->getFe()), *(advection->getQuadrature()),
			cellUpdateFlags);
	double result = 0.0;
	double area = 0.0;
	size_t dofs_per_cell = advection->getFe()->dofs_per_cell;

	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {

			// get global degrees of freedom
			std::vector<dealii::types::global_dof_index> localDoFIndices(
					dofs_per_cell);
			cell->get_dof_indices(localDoFIndices);
			// calculate the fe values for the cell
			feCellValues.reinit(cell);

			for (size_t i = 0; i < dofs_per_cell; i++) {
				result += ux(localDoFIndices.at(i))
						* feCellValues.JxW(
								advection->getCelldofToQIndex().at(i));
				area += feCellValues.JxW(advection->getCelldofToQIndex().at(i));
			}
		} /* if is locally owned */
	} /* for all cells */
	// communicate among MPI processes
	dealii::Utilities::MPI::MinMaxAvg global_res = dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	dealii::Utilities::MPI::MinMaxAvg global_area = dealii::Utilities::MPI::min_max_avg(area, MPI_COMM_WORLD);
	return global_res.sum / global_area.sum;

}
template double PhysicalProperties<2>::meanVelocityX(
		const distributed_vector& ux,
		shared_ptr<AdvectionOperator<2> > advection);
template double PhysicalProperties<3>::meanVelocityX(
		const distributed_vector& ux,
		shared_ptr<AdvectionOperator<3> > advection);

} /* namespace natrium */
