/*
 * ExponentialFilter.h
 *
 *  Created on: 01.02.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_SMOOTHING_EXPONENTIALFILTER_H_
#define LIBRARY_NATRIUM_SMOOTHING_EXPONENTIALFILTER_H_

#include "deal.II/fe/fe_dgq.h"
#include "deal.II/fe/fe_base.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/dofs/dof_handler.h"

#include "../utilities/BasicNames.h"
#include "Filter.h"

namespace natrium {

template <size_t dim>
class ExponentialFilter: public Filter<dim> {
private:
	double m_alpha;
	double m_s;

	size_t m_p;
	const dealii::Quadrature<dim>& m_quadrature;
	const dealii::FE_DGQ<dim>& m_sourceFE;

	// space of legendre polynomials (not using deal.II's FE classes, because the legendre fe
	// is inherited of fe_dgp, not fe_dqq. FE_DGP defines a COMPLETE polynomial space, meaning
	// that the basis is constructed by its maximum degree degree_x+degree_y+degree_z < p instead
	// of a tensor product with degree_x < p, degree_y <p, degree_z < p.

	std::vector<dealii::Polynomials::Polynomial<double> > m_legendre1D;

	numeric_matrix m_projectToLegendre;
	numeric_matrix m_projectFromLegendre;

	vector<size_t> m_degreeSum;
	vector<size_t> m_degreeMax;


	void makeProjectionMatrices(numeric_matrix& to_legendre, numeric_matrix& from_legendre);
	void makeDegreeVectors(size_t p);

	/**
	 * @short
	 *
	 */
	double evaluateLegendreND(size_t i, const dealii::Point<dim>& x);


public:
	ExponentialFilter(double alpha, double s, const dealii::Quadrature<dim>& quadrature, const dealii::FE_DGQ<dim>& fe);
	virtual ~ExponentialFilter(){};
	virtual void applyFilter(const dealii::DoFHandler<dim>& dof_handler, distributed_vector& dof_vector);

	double getAlpha() const {
		return m_alpha;
	}

	void setAlpha(double alpha) {
		m_alpha = alpha;
	}

	const numeric_matrix& getProjectFromLegendre() const {
		return m_projectFromLegendre;
	}

	const numeric_matrix& getProjectToLegendre() const {
		return m_projectToLegendre;
	}

	double getS() const {
		return m_s;
	}

	void setS(double s) {
		m_s = s;
	}

	const std::vector<dealii::Polynomials::Polynomial<double> >& getLegendre1D() const {
		return m_legendre1D;
	}
};

} /* namespace natrium */


#endif /* LIBRARY_NATRIUM_SMOOTHING_EXPONENTIALFILTER_H_ */
