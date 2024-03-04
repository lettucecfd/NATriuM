/*
 * ContinuousBoundaryGradient.cpp
 *
 *  Created on: 28.09.2016
 *      Author: akraem3m
 */

#include "ContinuousBoundaryGradient.h"

#ifndef NO_SUPPORT_HP_PARALLEL
#include "deal.II/hp/q_collection.h"
#include "deal.II/hp/mapping_collection.h"
#endif
#include "deal.II/numerics/vector_tools.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/numerics/fe_field_function.h"
#include "deal.II/fe/fe_nothing.h"
#include <deal.II/lac/affine_constraints.h>

namespace natrium {

template<size_t dim>
ContinuousBoundaryGradient<dim>::ContinuousBoundaryGradient(
		const dealii::DoFHandler<dim>& dof, const dealii::Mapping<dim>& mapping,
		const dealii::Quadrature<1>& quadrature,
		const dealii::Quadrature<dim - 1>& face_quadrature,
		const dealii::Quadrature<dim>& cell_quadrature,
		std::set<size_t> boundary_ids) :
		m_dof(dof), m_mapping(mapping), m_quadrature1D(quadrature), m_faceQuadrature(
				face_quadrature), m_cellQuadrature(cell_quadrature), m_boundaryIds(
				boundary_ids)
#ifdef NO_SUPPORT_HP_PARALLEL
				, m_FeDg(quadrature), m_FeConti(quadrature)
#endif
{
#ifndef NO_SUPPORT_HP_PARALLEL
	m_FeDgCollection.push_back(dealii::FESystem<dim>(dealii::FE_Nothing<dim>(),dim));
	dealii::FE_DGQArbitraryNodes<dim> fe_dgq(m_quadrature1D);
	m_FeDgCollection.push_back(dealii::FESystem<dim>(fe_dgq, dim));
	m_FeContiCollection.push_back(dealii::FESystem<dim>(dealii::FE_Nothing<dim>(),dim));
	m_FeContiCollection.push_back(dealii::FESystem<dim>(m_dof.get_fe(), dim));
#endif
}

template<size_t dim>
ContinuousBoundaryGradient<dim>::~ContinuousBoundaryGradient() {
	m_discontinuousBoundaryDoF->clear();
	m_continuousBoundaryDoF->clear();
}

template<size_t dim>
void ContinuousBoundaryGradient<dim>::reinit() {

	// set dofs for discontinuous boundaries
	m_discontinuousBoundaryDoF = boost::make_shared<DoFHandler<dim> >(
			m_dof.get_triangulation());

#ifndef NO_SUPPORT_HP_PARALLEL

	size_t faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
	for (typename DoFHandler<dim>::active_cell_iterator cell =
			m_discontinuousBoundaryDoF->begin_active();
			cell != m_discontinuousBoundaryDoF->end(); ++cell) {
		if (not cell->is_locally_owned()) {
			continue;
		}
		bool found = false;
		if (cell->at_boundary()) {
			for (size_t f = 0; f < faces_per_cell; f++) {
				if (m_boundaryIds.find(cell->face(f)->boundary_id())
						!= m_boundaryIds.end()) {
					cell->set_active_fe_index(1);
					found = true;
					continue;
				}
			}
		}
		if (not found) {
			cell->set_active_fe_index(0);
		}
	}
	m_discontinuousBoundaryDoF->distribute_dofs(m_FeDgCollection);
#else
	m_discontinuousBoundaryDoF->distribute_dofs(m_FeDg);
#endif

	m_continuousBoundaryDoF = boost::make_shared<DoFHandler<dim> >(
			m_dof.get_triangulation());
#ifndef NO_SUPPORT_HP_PARALLEL
	// set dofs for continuous boundaries
	for (typename DoFHandler<dim>::active_cell_iterator cell =
			m_continuousBoundaryDoF->begin_active();
			cell != m_continuousBoundaryDoF->end(); ++cell) {
		if (not cell->is_locally_owned()) {
			continue;
		}
		bool found = false;
		if (cell->at_boundary()) {
			for (size_t f = 0; f < faces_per_cell; f++) {
				if (m_boundaryIds.find(cell->face(f)->boundary_id())
						!= m_boundaryIds.end()) {
					cell->set_active_fe_index(1);
					found = true;
					continue;
				}
			}
		}
		if (not found) {
			cell->set_active_fe_index(0);
		}
	}
	m_continuousBoundaryDoF->distribute_dofs(m_FeContiCollection);
#else
	m_continuousBoundaryDoF->distribute_dofs(m_FeConti);
#endif

	m_discontinuousGradient.reinit(
			m_discontinuousBoundaryDoF->locally_owned_dofs());
	m_continuousGradient.reinit(m_continuousBoundaryDoF->locally_owned_dofs());
}

template<size_t dim>
void ContinuousBoundaryGradient<dim>::calculateGradients(
		const dealii::TrilinosWrappers::MPI::Vector& solution) {

	size_t faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
	dealii::UpdateFlags original_flags = dealii::update_gradients;
	dealii::UpdateFlags dg_flags = dealii::update_values;
	dealii::FEFaceValues<dim> orig_face_values(m_mapping, m_dof.get_fe(),
			m_faceQuadrature, original_flags);
#ifdef NO_SUPPORT_HP_PARALLEL
	FEFaceValues<dim> dg_face_values(m_mapping,
					m_discontinuousBoundaryDoF->get_fe(), m_faceQuadrature,
					dg_flags);
	size_t dg_dofs_per_cell =
			m_discontinuousBoundaryDoF->get_fe().dofs_per_cell;
#else
	dealii::hp::MappingCollection<dim> mapping_collection(m_mapping);
	dealii::hp::QCollection<dim - 1> face_quadrature_collection(
			m_faceQuadrature);
	dealii::hp::QCollection<dim> cell_quadrature_collection(m_cellQuadrature);
	dealii::hp::FEFaceValues<dim> dg_face_values(mapping_collection,
			m_discontinuousBoundaryDoF->get_fe(), face_quadrature_collection,
			dg_flags);
	size_t dg_dofs_per_cell =
			m_discontinuousBoundaryDoF->get_fe()[1].dofs_per_cell;
#endif
	std::vector<dealii::types::global_dof_index> dg_indices;
	dg_indices.resize(dg_dofs_per_cell);

	std::vector<dealii::Tensor<1, dim> > gradients;
	size_t n_q_points = m_faceQuadrature.size();
	gradients.resize(n_q_points);
	typename dealii::DoFHandler<dim>::active_cell_iterator original_cell =
			m_dof.begin_active();
	typename DoFHandler<dim>::active_cell_iterator dg_cell =
			m_discontinuousBoundaryDoF->begin_active();

	// iterate over both dof handlers, cells will come in the same order
	for (; original_cell != m_dof.end(); ++original_cell, ++dg_cell) {
		if (not original_cell->is_locally_owned()) {
			continue;
		}
		assert(dg_cell->n_active_fe_indices() == 1);
		if (dg_cell->active_fe_index() == 0) {
			continue;
		}
		for (size_t f = 0; f < faces_per_cell; f++) {
			if (m_boundaryIds.find(original_cell->face(f)->boundary_id())
					!= m_boundaryIds.end()) {
				continue;
			}
			// calculate gradients
			dg_cell->get_dof_indices(dg_indices);
			orig_face_values.reinit(original_cell, f);
			dg_face_values.reinit(dg_cell, f);
			orig_face_values.get_function_gradients(solution, gradients);

			for (size_t i = 0; i < dg_dofs_per_cell; i++) {

#ifdef NO_SUPPORT_HP_PARALLEL
				size_t component_i =
						m_discontinuousBoundaryDoF->get_fe().system_to_component_index(
								i).first;
				size_t index_i =
						m_discontinuousBoundaryDoF->get_fe().system_to_component_index(
								i).second;
#else
				size_t component_i =
						m_discontinuousBoundaryDoF->get_fe()[dg_cell->active_fe_index()].system_to_component_index(
								i).first;
				size_t index_i =
						m_discontinuousBoundaryDoF->get_fe()[dg_cell->active_fe_index()].system_to_component_index(
								i).second;
#endif
				size_t q_point = 0;
				// select the right quadrature point for this dof
				for (; q_point < n_q_points; q_point++) {
					double shape =
							dg_face_values.get_present_fe_values().shape_value_component(
									index_i, q_point, component_i);
					if (shape > 1e-20) {
						// assert shape == 1
						assert(1 + 1e-20 > shape);
						assert(1 - 1e-20 < shape);
						break;
					}
				}
				m_discontinuousGradient(dg_indices.at(i)) = gradients.at(
						index_i)[component_i];
			}
		}
	}

	// project from discontinuous to continuous fe
#ifdef NO_SUPPORT_HP_PARALLEL
//TODO BRUTE FORCE KILLED SEDG BOUNDARY
/*dealii::Functions::FEFieldFunction<dim> fe_function(
			*m_discontinuousBoundaryDoF, m_discontinuousGradient, m_mapping);
*/
#else
	dealii::Functions::FEFieldFunction<dim, typename DoFHandler<dim> > fe_function(
			*m_discontinuousBoundaryDoF, m_discontinuousGradient, m_mapping);
#endif
	dealii::AffineConstraints<double> c;
	c.close();
#ifdef NO_SUPPORT_HP_PARALLEL
/*	dealii::VectorTools::project(m_mapping, *m_continuousBoundaryDoF,
			c, m_cellQuadrature, fe_function, m_continuousGradient,
			false, m_faceQuadrature); */
	m_feFaceValues = boost::make_shared<FEFaceValues<dim> >(
					m_mapping, m_continuousBoundaryDoF->get_fe(),
					m_faceQuadrature, dg_flags);
#else
	dealii::VectorTools::project(mapping_collection, *m_continuousBoundaryDoF,
			c, cell_quadrature_collection, fe_function, m_continuousGradient,
			false, face_quadrature_collection);	// prepare face values object for accessing gradients
	m_feFaceValues = boost::make_shared<dealii::hp::FEFaceValues<dim> >(
			mapping_collection, m_continuousBoundaryDoF->get_fe(),
			face_quadrature_collection, dg_flags);
#endif


	m_dofIndizes.clear();
	m_dofIndizes.resize(m_continuousBoundaryDoF->get_fe(1).dofs_per_cell);
}

template<size_t dim>
double ContinuousBoundaryGradient<dim>::get_gradient_component(size_t q_point,
		size_t component) {
	size_t dofs_per_cell = m_dofIndizes.size();
	for (size_t i = 0; i < dofs_per_cell; i++) {
#ifdef NO_SUPPORT_HP_PARALLEL
		size_t component_i =
				m_continuousBoundaryDoF->get_fe().system_to_component_index(
						i).first;
#else
		size_t component_i =
				m_continuousBoundaryDoF->get_fe()[1].system_to_component_index(
						i).first;
#endif
		if (component != component_i) {
			continue;
		}
#ifdef NO_SUPPORT_HP_PARALLEL
		size_t index_i =
				m_continuousBoundaryDoF->get_fe().system_to_component_index(
						i).second;
#else
		size_t index_i =
				m_continuousBoundaryDoF->get_fe()[1].system_to_component_index(
						i).second;
#endif
		double shape =
				m_feFaceValues->get_present_fe_values().shape_value_component(
						index_i, q_point, component_i);
		if (shape > 1e-20) {
			// assert shape == 1
			assert(1 + 1e-20 > shape);
			assert(1 - 1e-20 < shape);

			return m_continuousGradient(m_dofIndizes.at(i));
		}
	}
	// not found: throw exception
	assert(false);
	return 0;
}
template class ContinuousBoundaryGradient<2> ;
template class ContinuousBoundaryGradient<3> ;

} /* namespace natrium */
