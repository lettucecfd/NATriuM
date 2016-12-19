/*
 * AdvectionTools.h
 *
 *  Created on: 26.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_ADVECTION_ADVECTIONTOOLS_H_
#define LIBRARY_NATRIUM_ADVECTION_ADVECTIONTOOLS_H_

#include "deal.II/base/quadrature.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/fe/fe.h"
#include "../solver/SolverConfiguration.h"
#include "../utilities/BasicNames.h"


namespace natrium {

namespace AdvectionTools{

/**
 * @short create a quadrature instance by its name
 * @param name The quadrature name (as an enum, defined in ConfigNames)
 * @param n_points The number of points in the quadrature (in one dimension)
 * @return a shared pointer to the quadrature
 */
template <size_t dim>
boost::shared_ptr<dealii::Quadrature<dim> > make_quadrature_by_name(QuadratureName name, size_t n_points);

/**
 * @short create a Lagrangian finite element instance by its name
 * @param name The finite element name (as an enum, defined in ConfigNames)
 * @param order the order of the 1D polynomial shape functions
 * @param dg if true, this function returns a discontinuous Galerkin-type finite element (FE_DGQArbitraryNodes),
 * 				else, it returns a continuous Galerkin-type finite element (FE_Q)
 * @return a shared pointer to the finite element
 */
template <size_t dim>
boost::shared_ptr<dealii::FiniteElement<dim> > make_fe_by_name(SupportPointsName name, size_t order, bool dg);

} /* namespace AdvectionTools */


/**
 * @short Exception class for AdvectionOperator
 */
class AdvectionSolverException: public NATriuMException {
private:
	std::string message;
public:
	AdvectionSolverException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	AdvectionSolverException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~AdvectionSolverException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_ADVECTIONTOOLS_H_ */
