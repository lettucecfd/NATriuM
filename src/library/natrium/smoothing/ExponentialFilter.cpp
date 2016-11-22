/*
 * ExponentialFilter.cpp
 *
 *  Created on: 01.02.2016
 *      Author: akraem3m
 */

#include "ExponentialFilter.h"
#include "deal.II/fe/fe_tools.h"
#include "deal.II/base/polynomial_space.h"
#include "deal.II/base/polynomial.h"

namespace natrium {

template<size_t dim>
ExponentialFilter<dim>::ExponentialFilter(double alpha, double s, size_t Nc,
		bool by_sum, const dealii::Quadrature<dim>& quadrature,
		const dealii::FiniteElement<dim>& fe) :
		m_alpha(alpha), m_s(s), m_Nc(Nc), m_bySum(by_sum), m_p(fe.degree), m_quadrature(
				quadrature), m_sourceFE(fe), m_legendre1D(0), m_projectToLegendre(
				pow(m_sourceFE.degree + 1, dim)), m_projectFromLegendre(
				pow(m_sourceFE.degree + 1, dim)) {
	m_legendre1D = dealii::Polynomials::Legendre::generate_complete_basis(m_p);
	makeProjectionMatrices(m_projectToLegendre, m_projectFromLegendre);
	makeDegreeVectors(m_p);

}

template<size_t dim>
void natrium::ExponentialFilter<dim>::makeProjectionMatrices(
		numeric_matrix& to_legendre, numeric_matrix& from_legendre) {

	// assert right sizes
	const size_t n = to_legendre.n();
	assert(n != 0);
	assert(n == to_legendre.m());
	assert(n == from_legendre.n());
	assert(n == from_legendre.m());
	assert(n == pow(m_p + 1, dim));

	// fill matrices
	const std::vector<dealii::Point<dim> > & quad_points =
			m_quadrature.get_points();
	const std::vector<double> & weights = m_quadrature.get_weights();
	size_t Q = quad_points.size();
	assert(Q == weights.size());
	dealii::FullMatrix<double> quad_source(n);
	dealii::FullMatrix<double> quad_legendre(n);

	for (size_t k = 0; k < n; k++) {
		for (size_t i = 0; i < n; i++) {
			for (size_t q = 0; q < Q; q++) {
				double phi_i = m_sourceFE.shape_value(i, quad_points.at(q));
				double psi_i = evaluateLegendreND(i, quad_points.at(q));
				double psi_k = evaluateLegendreND(k, quad_points.at(q));
				quad_source(k, i) = quad_source(k, i)
						+ (weights.at(q) * phi_i * psi_k);
				quad_legendre(k, i) = quad_legendre(k, i)
						+ (weights.at(q) * psi_i * psi_k);
			}
		}
	}
	dealii::FullMatrix<double> quad_legendre_inv(n);
	quad_legendre_inv.invert(quad_legendre);
	quad_legendre_inv.mmult(m_projectToLegendre, quad_source);
	m_projectFromLegendre.invert(m_projectToLegendre);

}

template<>
inline double natrium::ExponentialFilter<1>::evaluateLegendreND(size_t i,
		const dealii::Point<1>& x) {
	assert(i < m_legendre1D.size());
	return m_legendre1D.at(i).value(x(0));
}
template<>
inline double natrium::ExponentialFilter<2>::evaluateLegendreND(size_t i,
		const dealii::Point<2>& x) {
	assert(i < (m_p + 1) * (m_p + 1));
	size_t ix = i / (m_p + 1);
	size_t iy = i % (m_p + 1);
	return m_legendre1D.at(ix).value(x(0)) * m_legendre1D.at(iy).value(x(1));
}
template<>
inline double natrium::ExponentialFilter<3>::evaluateLegendreND(size_t i,
		const dealii::Point<3>& x) {
	assert(i < (m_p + 1) * (m_p + 1) * (m_p + 1));
	size_t ix = i / ((m_p + 1) * (m_p + 1));
	size_t iy = (i % ((m_p + 1) * (m_p + 1))) / (m_p + 1);
	size_t iz = i % (m_p + 1);
	return m_legendre1D.at(ix).value(x(0)) * m_legendre1D.at(iy).value(x(1))
			* m_legendre1D.at(iz).value(x(2));
}

template<>
inline void ExponentialFilter<1>::makeDegreeVectors(size_t p) {
	m_degreeMax.resize(p + 1);
	m_degreeSum.resize(p + 1);
	for (size_t i = 0; i < p + 1; i++) {
		m_degreeMax.at(i) = i;
		m_degreeSum.at(i) = i;
	}
}

template<>
inline void ExponentialFilter<2>::makeDegreeVectors(size_t p) {
	m_degreeMax.resize((p + 1) * (p + 1));
	m_degreeSum.resize((p + 1) * (p + 1));
	for (size_t i = 0; i < (p + 1) * (p + 1); i++) {
		size_t ix = i / (p + 1);
		size_t iy = i % (p + 1);
		size_t max = ix;
		if (iy > max)
			max = iy;
		m_degreeMax.at(i) = max;
		m_degreeSum.at(i) = ix + iy;
	}
}

template<>
inline void ExponentialFilter<3>::makeDegreeVectors(size_t p) {
	m_degreeMax.resize((p + 1) * (p + 1) * (p + 1));
	m_degreeSum.resize((p + 1) * (p + 1) * (p + 1));
	for (size_t i = 0; i < (p + 1) * (p + 1) * (p + 1); i++) {
		size_t ix = i / ((p + 1) * (p + 1));
		size_t max = ix;
		size_t iy = ( i % (p+1) * (p+1) ) / (p + 1);
		if (iy > max)
			max = iy;
		size_t iz = i % (p + 1);
		if (iz > max)
			max = iz;
		m_degreeMax.at(i) = max;
		m_degreeSum.at(i) = ix + iy + iz;
	}
}

template<size_t dim>
void ExponentialFilter<dim>::applyFilter(
		const dealii::DoFHandler<dim>& dof_handler,
		distributed_vector& dof_vector) {

	const size_t dofs_per_cell = m_sourceFE.dofs_per_cell;
	std::vector<dealii::types::global_dof_index> local_dof_indices(
			dofs_per_cell);
	numeric_vector source_fe_dofs(dofs_per_cell);
	numeric_vector legendre_dofs(dofs_per_cell);
	size_t max_degree = m_p;
	if (m_bySum)
		max_degree = dim * m_p;
	double sigma = 0.0; // filter scaling

	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			// get global degrees of freedom
			cell->get_dof_indices(local_dof_indices);

			// project to legendre basis
			for (size_t i = 0; i < dofs_per_cell; i++) {
				source_fe_dofs(i) = dof_vector(local_dof_indices[i]);
			}
			m_projectToLegendre.vmult(legendre_dofs, source_fe_dofs);

			// apply exponential filtering
			for (size_t i = 0; i < dofs_per_cell; i++) {
				size_t degree = m_degreeMax.at(i);
				if (m_bySum)
					degree = m_degreeSum.at(i);

				if (degree >= m_Nc) {
					sigma = std::exp(
							-m_alpha
									* pow(
											((double) degree + 1 - m_Nc)
													/ ((double) max_degree + 1
															- m_Nc), m_s));
					legendre_dofs(i) = sigma * legendre_dofs(i);
				}
			}

			// project back to original basis
			m_projectFromLegendre.vmult(source_fe_dofs, legendre_dofs);
			for (size_t i = 0; i < dofs_per_cell; i++) {
				dof_vector(local_dof_indices[i]) = source_fe_dofs(i);
			}

		} /* if is locally owned */

	} /* for cells */

} /* applyFilter */

template class ExponentialFilter<1> ;
template class ExponentialFilter<2> ;
template class ExponentialFilter<3> ;

} /* namespace natrium */
