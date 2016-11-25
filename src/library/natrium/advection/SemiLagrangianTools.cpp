/*
 * SemiLagrangianTools.cpp
 *
 *  Created on: 17.11.2016
 *      Author: akraem3m
 */

#include <array>

#include "SemiLagrangianTools.h"

namespace natrium {

template<size_t dim>
void shapeFunctionValue(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const std::vector<dealii::Point<dim> >& points,
		std::vector<std::vector<double> >&values,
		const dealii::Mapping<dim>& mapping) {

	//TimerOutput::Scope timer_section(Timing::getTimer(), "Assembly: evaluate function");

	const size_t n_dofs_per_cell = cell->get_fe().n_dofs_per_cell();
	const size_t n_points = points.size();
	assert(n_points > 0);
	assert(values.size() == n_points);

	// clear values
	for (size_t i = 0; i < n_points; i++) {
		assert(values.at(i).size() == n_dofs_per_cell);
		for (size_t j = 0; j < n_dofs_per_cell; j++) {
			values[i][j] = 0;
		}
	}
	// transform to unit points
	std::vector<dealii::Point<dim> > unit_points;
	//std::vector<double> weights (n_points, 1.0);
	for (size_t i = 0; i < n_points; i++) {
		dealii::Point<dim> h = mapping.transform_real_to_unit_cell(cell,
				points[i]);
		for (size_t i = 0; i < dim; i++) {
			if (fabs(h[i]) < 1e-10) {
				h[i] = 0;
			}
			if (fabs(h[i] - 1) < 1e-10) {
				h[i] = 1;
			}
			assert(h[i] <= 1);
			assert(h[i] >= 0);
		}
		unit_points.push_back(h);
	}
	// Now we can find out about the points
	for (size_t i = 0; i < n_points; i++) {
		dealii::Quadrature<dim> quad(unit_points.at(i)); //, weights);
		dealii::FEValues<dim> fe_v(mapping, cell->get_fe(), quad,
				dealii::update_values);
		fe_v.reinit(cell);
		for (size_t j = 0; j < n_dofs_per_cell; j++) {
			values[i][j] = fe_v.shape_value(j, 0);
		}
	}
}

// explicit template instantiation
template void shapeFunctionValue<2>(
		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		const std::vector<dealii::Point<2> >& points,
		std::vector<std::vector<double> >&values,
		const dealii::Mapping<2>& mapping);
template void shapeFunctionValue<3>(
		const typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		const std::vector<dealii::Point<3> >& points,
		std::vector<std::vector<double> >&values,
		const dealii::Mapping<3>& mapping);

template<size_t dim>
int supportPointNr(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const dealii::Point<dim>& point, const dealii::Quadrature<dim>& quad,
		const dealii::MappingQ<dim>& mapping) {

	// Now we can find out about the points
	dealii::FEValues<dim> fe_v(mapping, cell->get_fe(), quad,
			dealii::update_quadrature_points);
	fe_v.reinit(cell);
	for (size_t q = 0; q < quad.size(); q++) {
		if (fe_v.quadrature_point(q).distance(point) < 1e-10){
			return q;
		}
	}
	return -1;
}
template
int supportPointNr<2>(
		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		const dealii::Point<2>& point, const dealii::Quadrature<2>& quad,
		const dealii::MappingQ<2>& mapping);
template
int supportPointNr<3>(
		const typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		const dealii::Point<3>& point, const dealii::Quadrature<3>& quad,
		const dealii::MappingQ<3>& mapping);

template<size_t dim>
boost::shared_ptr<dealii::FEValues<dim> > reinitArbitraryPoints(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const std::vector<dealii::Point<dim> >& points,
		const dealii::Mapping<dim>& mapping, const dealii::UpdateFlags& flags) {

	//TimerOutput::Scope timer_section(Timing::getTimer(), "Assembly: evaluate function");
	const size_t n_points = points.size();
	assert(n_points > 0);

	// transform to unit points
	std::vector<dealii::Point<dim> > unit_points;
	//std::vector<double> weights (n_points, 1.0);
	for (size_t i = 0; i < n_points; i++) {
		dealii::Point<dim> h = mapping.transform_real_to_unit_cell(cell,
				points[i]);
		for (size_t i = 0; i < dim; i++) {
			if (fabs(h[i]) < 1e-10) {
				h[i] = 0;
			}
			if (fabs(h[i] - 1) < 1e-10) {
				h[i] = 1;
			}
			assert(h[i] <= 1);
			assert(h[i] >= 0);
		}
		unit_points.push_back(h);
	}
	// Now we can find out about the points
	dealii::Quadrature<dim> quad(unit_points); //, weights);
	boost::shared_ptr<dealii::FEValues<dim> > fe_v = boost::make_shared<
			dealii::FEValues<dim> >(mapping, cell->get_fe(), quad, flags);
	fe_v->reinit(cell);
	return fe_v;
}

template boost::shared_ptr<dealii::FEValues<2> > reinitArbitraryPoints<2>(
		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		const std::vector<dealii::Point<2> >& points,
		const dealii::Mapping<2>& mapping, const dealii::UpdateFlags& flags);
template boost::shared_ptr<dealii::FEValues<3> > reinitArbitraryPoints<3>(
		const typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		const std::vector<dealii::Point<3> >& points,
		const dealii::Mapping<3>& mapping, const dealii::UpdateFlags& flags);




} /* namespace natrium */
