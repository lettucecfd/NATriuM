/*
 * AdvectionTools.cpp
 *
 *  Created on: 26.11.2016
 *      Author: akraem3m
 */

#include "AdvectionTools.h"
#include "../utilities/NATriuMException.h"

#include "deal.II/fe/fe_dgq.h"
#include "deal.II/fe/fe_q.h"

namespace natrium {

namespace AdvectionTools {

template<size_t dim>
boost::shared_ptr<dealii::Quadrature<dim> > make_quadrature_by_name(
		QuadratureName name, size_t n_points) {
	switch (name) {
	case QGAUSS_LOBATTO: {
		return boost::make_shared<dealii::QGaussLobatto<dim> >(n_points);
		break;
	}
	case QGAUSS: {
		return boost::make_shared<dealii::QGauss<dim> >(n_points);
		break;
	}
	default: {
		throw NATriuMException(
				"The quadrature you specified in your solver configuration"
						"was not supported by AdvectionTools::make_quadrature_by_name()");
		break;
	}
	}
}

// explicit templates
template boost::shared_ptr<dealii::Quadrature<1> > make_quadrature_by_name<1>(
		QuadratureName name, size_t n_points);
template boost::shared_ptr<dealii::Quadrature<2> > make_quadrature_by_name<2>(
		QuadratureName name, size_t n_points);
template boost::shared_ptr<dealii::Quadrature<3> > make_quadrature_by_name<3>(
		QuadratureName name, size_t n_points);

template<size_t dim>
boost::shared_ptr<dealii::FiniteElement<dim> > make_fe_by_name(
		SupportPointsName name, size_t order, bool dg) {

	// make base quadrature
	boost::shared_ptr<dealii::Quadrature<1> > quad;
	switch (name) {
	case GAUSS_LOBATTO_POINTS: {
		quad = boost::make_shared<dealii::QGaussLobatto<1> >(order + 1);
		break;
	}
	case GAUSS_CHEBYSHEV_POINTS: {
		quad = boost::make_shared<dealii::QGaussChebyshev<1> >(order + 1);
		break;
	}
	case GAUSS_LOBATTO_CHEBYSHEV_POINTS: {
		quad = boost::make_shared<dealii::QGaussLobattoChebyshev<1> >(order + 1);
		break;
	}
	case EQUIDISTANT_POINTS: {
		quad = boost::make_shared<dealii::QGauss<1> >(order + 1);
		break;
	}
	default: {
		throw NATriuMException(
				"The set of support points that you specified in your solver configuration"
						"was not supported by AdvectionTools::make_fe_by_name()");
		break;
	}
	}

	// make fe
	if (dg){
		return boost::make_shared<dealii::FE_DGQArbitraryNodes<dim> >(*quad);
		if (name != GAUSS_LOBATTO_POINTS){
			LOG(WARNING) << "Using another FE than GaussLobatto in connection with the DG method is "
					"highly discouraged as it would require to solve a linear equation system "
					"in each time step. This is very costly and thus not implemented in NATriuM."
					<< endl;
		}
	} else {
		return boost::make_shared<dealii::FE_Q<dim> >(*quad);
	}
}

template
boost::shared_ptr<dealii::FiniteElement<1> > make_fe_by_name<1>(
		SupportPointsName name, size_t order, bool dg);
template
boost::shared_ptr<dealii::FiniteElement<2> > make_fe_by_name<2>(
		SupportPointsName name, size_t order, bool dg);
template
boost::shared_ptr<dealii::FiniteElement<3> > make_fe_by_name<3>(
		SupportPointsName name, size_t order, bool dg);

} /* namespace AdvectionTools */

} /* namespace natrium */
