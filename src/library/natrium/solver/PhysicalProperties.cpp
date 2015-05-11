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
		const vector<distributed_vector>& u, const distributed_vector& rho) {
	const size_t n_dofs = u.at(0).size();
	assert(n_dofs == rho.size());
	assert(n_dofs == u.at(1).size());
	double kE = 0.0;
	for (size_t i = 0; i < n_dofs; i++) {
		kE += 0.5 * rho(i)
				* (u.at(0)(i) * u.at(0)(i) + u.at(1)(i) * u.at(1)(i));
	}
	return kE * 0.5 / n_dofs;
}
template<> double PhysicalProperties<3>::kineticEnergy(
		const vector<distributed_vector>& u, const distributed_vector& rho) {
	const size_t n_dofs = u.at(0).size();
	assert(n_dofs == rho.size());
	assert(n_dofs == u.at(1).size());
	assert(n_dofs == u.at(2).size());
	double kE = 0.0;
	for (size_t i = 0; i < n_dofs; i++) {
		kE += rho(i)
				* (u.at(0)(i) * u.at(0)(i) + u.at(1)(i) * u.at(1)(i)
						+ u.at(2)(i) * u.at(2)(i));
	}
	return kE * 0.5 / n_dofs;
}

template<size_t dim>
double PhysicalProperties<dim>::maximalPressure(const distributed_vector& rho,
		const double speedOfSound, double & minimalPressure) {
	const size_t n_dofs = rho.size();
	double maximal = -100000000000.;
	minimalPressure = 10000000000000000;
	for (size_t i = 0; i < n_dofs; i++) {
		if (rho(i) < minimalPressure) {
			minimalPressure = rho(i);
		}
		if (rho(i) > maximal) {
			maximal = rho(i);
		}
	}
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
double PhysicalProperties<dim>::massFluxX(const distributed_vector& ux,
		shared_ptr<AdvectionOperator<dim> > advection, double Lx) {

	// Integrate ux over whole domain
	const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values;
	const dealii::DoFHandler<dim> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<dim> feCellValues(advection->getMapping(), *(advection->getFe()), *(advection->getQuadrature()),
				cellUpdateFlags);
	double result = 0.0;
	size_t dofs_per_cell = advection->getFe()->dofs_per_cell;


	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {

		// get global degrees of freedom
		std::vector<dealii::types::global_dof_index> localDoFIndices(dofs_per_cell);
		cell->get_dof_indices(localDoFIndices);
		// calculate the fe values for the cell
		feCellValues.reinit(cell);

		for (size_t i = 0; i < dofs_per_cell; i++){
			result += ux(localDoFIndices.at(i)) * feCellValues.JxW(advection->getCelldofToQIndex().at(i));
		}
	}
	return result / Lx;
}
template double PhysicalProperties<2>::massFluxX(const distributed_vector& ux,
		shared_ptr<AdvectionOperator<2> > advection, double Lx);
template double PhysicalProperties<3>::massFluxX(const distributed_vector& ux,
		shared_ptr<AdvectionOperator<3> > advection, double Lx);

} /* namespace natrium */
