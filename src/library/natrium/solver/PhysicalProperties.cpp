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
		boost::shared_ptr<AdvectionOperator<2> > advection) {
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
	dealii::Utilities::MPI::MinMaxAvg global_res =
			dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	return 0.5 * global_res.sum;

}
template<> double PhysicalProperties<3>::kineticEnergy(
		const vector<distributed_vector>& u, const distributed_vector& rho,
		boost::shared_ptr<AdvectionOperator<3> > advection) {
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
	dealii::Utilities::MPI::MinMaxAvg global_res =
			dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	return 0.5 * global_res.sum;
}

template<> double PhysicalProperties<2>::enstrophy(
		const vector<distributed_vector>& u,
		boost::shared_ptr<AdvectionOperator<2> > advection) {
	const size_t n_dofs = u.at(0).size();
	assert(n_dofs == u.at(1).size());

	const distributed_vector& ux = u.at(0);
	const distributed_vector& uy = u.at(1);

	// Integrate ux over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values
			| dealii::update_gradients;
	const dealii::DoFHandler<2> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<2> feCellValues(advection->getMapping(),
			*(advection->getFe()), *(advection->getQuadrature()),
			cellUpdateFlags);
	double result = 0.0;
	double vorticity = 0.0;
	size_t local_i;
	size_t q_point;
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

			// Enstrophy = int w^2 = int ( dv/dx - du/dy)^2  = w_q ( v_i dphi/dx (x_q) - u_j dphi/dy (x_q) )
			for (size_t q = 0; q < dofs_per_cell; q++) {
				vorticity = 0.0;
				q_point = advection->getCelldofToQIndex().at(q);
				for (size_t i = 0; i < dofs_per_cell; i++) {
					local_i = localDoFIndices.at(i);
					vorticity += (uy(local_i)
							* feCellValues.shape_grad(i, q_point)[0]
							- ux(local_i)
									* feCellValues.shape_grad(i, q_point)[1]);
				}
				result += vorticity * vorticity * feCellValues.JxW(q_point);
			}
		} /* if is locally owned */
	} /* for cells */

// communicate among MPI processes
	dealii::Utilities::MPI::MinMaxAvg global_res =
			dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	// the enstrophy can be defined with and without factor 1/2; here: without
	return global_res.sum;

}

template<> double PhysicalProperties<3>::enstrophy(
		const vector<distributed_vector>& u,
		boost::shared_ptr<AdvectionOperator<3> > advection) {
	const size_t n_dofs = u.at(0).size();
	assert(n_dofs == u.at(1).size());
	assert(n_dofs == u.at(2).size());

	const distributed_vector& ux = u.at(0);
	const distributed_vector& uy = u.at(1);
	const distributed_vector& uz = u.at(2);

	// Integrate ux over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values | dealii::update_gradients;
	const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<3> feCellValues(advection->getMapping(),
			*(advection->getFe()), *(advection->getQuadrature()),
			cellUpdateFlags);
	double result = 0.0;
	double frobenius_sq = 0.0;
	size_t local_i;
	size_t q_point;
	double du_dxk;
	double dv_dxk;
	double dw_dxk;
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

			// Enstrophy = int |frob(u)|^2 = int sum_k (ui dphi_i/dxk)^2 +  (vi dphi_i/dxk)^2 + (wi dphi_i/dxk)^2
			// = w_q  sum_k (ui dphi_i/dxk(x_q))^2 +  (vi dphi_i/dxk(x_q))^2 + (wi dphi_i/dxk(x_q))^2
			for (size_t q = 0; q < dofs_per_cell; q++) {
				frobenius_sq = 0.0;
				q_point = advection->getCelldofToQIndex().at(q);
				for (size_t k = 0; k < 3; k++) {
					du_dxk = 0;
					dv_dxk = 0;
					dw_dxk = 0;
					for (size_t i = 0; i < dofs_per_cell; i++) {
						local_i = localDoFIndices.at(i);
						du_dxk += ux(local_i)
								* feCellValues.shape_grad(i, q_point)[k];
						dv_dxk += uy(local_i)
								* feCellValues.shape_grad(i, q_point)[k];
						dw_dxk += uz(local_i)
								* feCellValues.shape_grad(i, q_point)[k];
					}
					frobenius_sq += du_dxk * du_dxk + dv_dxk * dv_dxk
							+ dw_dxk * dw_dxk;
				}
				result += frobenius_sq * feCellValues.JxW(q_point);
			}
		} /* if is locally owned */
	} /* for cells */

	// communicate among MPI processes
	dealii::Utilities::MPI::MinMaxAvg global_res =
			dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	return global_res.sum;
}

template<size_t dim>
double PhysicalProperties<dim>::mass(const distributed_vector& rho,
		boost::shared_ptr<AdvectionOperator<dim> > advection) {
	const size_t n_dofs = advection->getNumberOfDoFs();
	assert(n_dofs == rho.size());

	// Integrate rho over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values;
	const dealii::DoFHandler<dim> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<dim> feCellValues(advection->getMapping(),
			*(advection->getFe()), *(advection->getQuadrature()),
			cellUpdateFlags);
	double result = 0.0;
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

			size_t local_i;
			for (size_t i = 0; i < dofs_per_cell; i++) {
				local_i = localDoFIndices.at(i);
				result += rho(local_i)
						* feCellValues.JxW(
								advection->getCelldofToQIndex().at(i));
			}
		} /* if is locally owned */
	} /* for cells */

	// communicate among MPI processes
	dealii::Utilities::MPI::MinMaxAvg global_res =
			dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	return global_res.sum;
}

template<size_t dim>
double PhysicalProperties<dim>::maximalPressure(const distributed_vector& rho,
		const double speedOfSound, double & minimalPressure) {
	double maximal = rho.max();
	minimalPressure = rho.min();
	minimalPressure = (minimalPressure - 1.0) * speedOfSound * speedOfSound;
	return (maximal - 1.0) * (speedOfSound * speedOfSound);

}

template<size_t dim>
double PhysicalProperties<dim>::meanVelocityX(const distributed_vector& ux,
		boost::shared_ptr<AdvectionOperator<dim> > advection) {

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
	dealii::Utilities::MPI::MinMaxAvg global_res =
			dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
	dealii::Utilities::MPI::MinMaxAvg global_area =
			dealii::Utilities::MPI::min_max_avg(area, MPI_COMM_WORLD);
	return global_res.sum / global_area.sum;

}

// Explicit instantiation
template class PhysicalProperties<2>;
template class PhysicalProperties<3>;

} /* namespace natrium */
