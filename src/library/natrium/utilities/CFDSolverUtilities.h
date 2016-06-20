/**
 * @file CFDSolverUtilities.h
 * @short General utility functions for the CFD solver
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVERUTILITIES_H_
#define CFDSOLVERUTILITIES_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/quadrature_lib.h"

#include "../stencils/Stencil.h"
#include "../stencils/D2Q9.h"
#include "../stencils/D3Q15.h"
#include "../stencils/D3Q19.h"
#include "../stencils/D3Q27.h"

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "../solver/SolverConfiguration.h"

namespace natrium {

/** @short
 * Some tools for the CFDSolver and its simple usage.
 */
namespace CFDSolverUtilities {

class CFDSolverUtilitiesException: public NATriuMException {
private:
	std::string message;
public:
	CFDSolverUtilitiesException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	CFDSolverUtilitiesException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~CFDSolverUtilitiesException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short
 */
template<size_t dim>
double getMinimumDoFDistanceGLL(const Mesh<dim>& tria,
		const size_t orderOfFiniteElement);

template<size_t dim>
double getMinimumVertexDistance(const Mesh<dim>& tria);

template<size_t dim>
double calculateTimestep(const Mesh<dim>& tria,
		const size_t orderOfFiniteElement, const Stencil& stencil, double cFL =
				0.4);

/**
 * @short stolen from Deal.II's step 49 tutorial
 */
template<int dim>
void mesh_info(const Mesh<dim> &tria, const std::string &filename);

/**
 * @short Select an integrator by a specified id:
 *  - 1) RUNGE_KUTTA_5STAGE
 *  - 2) THETA_METHOD
 *  - 3) EXPONENTIAL
 *  - 4) FORWARD_EULER
 *  - 5) RK_THIRD_ORDER
 *  - 6) RK_CLASSIC_FOURTH_ORDER
 *  - 7) BACKWARD_EULER
 *  - 8) IMPLICIT_MIDPOINT
 *  - 9) CRANK_NICOLSON
 *  - 10) SDIRK_TWO_STAGES
 *  - 11) HEUN_EULER
 *  - 12) BOGACKI_SHAMPINE
 *  - 13) DOPRI
 *  - 14) FEHLBERG
 *  - 15) CASH_KARP
 *  @param[in] id Number of integrator
 *  @param[out] time_integrator Argument to be used in configuration->setTimeIntegrator()
 *  @param[out] deal_integrator Argument to be used in configuration->setDealIntegrator()
 *  @param[out] integrator_name The name of the integrator as a string.
 *  @note Integrators 4-15 are Deal.II's built-in integrators
 */
void get_integrator_by_id(size_t id, TimeIntegratorName& time_integrator,
		DealIntegratorName& deal_integrator, std::string& integrator_name);

boost::shared_ptr<Stencil> make_stencil(size_t d, size_t q, size_t scaling);


} /* CFDSolverUtilities */
} /* namespace natrium */

#endif /* CFDSOLVERUTILITIES_H_ */
