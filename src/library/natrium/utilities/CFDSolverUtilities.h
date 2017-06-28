/**
 * @file CFDSolverUtilities.h
 * @short General utility functions for the CFD solver
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVERUTILITIES_H_
#define CFDSOLVERUTILITIES_H_

#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe.h"
#include "deal.II/base/quadrature_lib.h"

#include "../stencils/Stencil.h"
#include "../stencils/D2Q9.h"
#include "../stencils/D3Q15.h"
#include "../stencils/D3Q19.h"
#include "../stencils/D3Q27.h"

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "../solver/SolverConfiguration.h"
#include "../utilities/Timing.h"

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
double getMinimumDoFDistance(const Mesh<dim>& tria,
		const dealii::FiniteElement<dim,dim>& fe);

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

std::string get_integrator_name(
		const TimeIntegratorName& time_integrator,
		const DealIntegratorName& deal_integrator) ;


boost::shared_ptr<Stencil> make_stencil(size_t d, size_t q, size_t scaling);



/**
 * @short get a writeable copy of the velocity
 *			The background of this function is that the velocity is stored in a ghosted vector (see deal.II glossary entry
 * 			on ghosted vector, when we use non-DG type discretizations). This is required, e.g. when integrating the
 * 			velocity on a cell, as all locally relevant dofs have to be stored on the local processor. Ghosted vectors
 * 			in deal.II can only be written to by assigning a non-ghosted vector. In order to write
 * 			something into the velocity vectors, we have to copy the values into non-ghosted vectors with the present
 * 			function, change then and assign the writeable vector to the member variable via applyWriteableVeloctiy()).
 *
 */
inline std::vector<distributed_vector>& getWriteableVelocity(std::vector<distributed_vector>& writeable, const std::vector<distributed_vector>& member, const dealii::IndexSet& locally_owned){

	TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
	writeable.resize(member.size());
	for (size_t i = 0; i < member.size(); i++){
		writeable.at(i).reinit(locally_owned, member.at(i).get_mpi_communicator(), true);
		assert (not writeable.at(i).has_ghost_elements());
		writeable.at(i) = member.at(i);
		assert (not writeable.at(i).has_ghost_elements());
	}
	return writeable;
}

/**
 * @short copy the changes in the writeable copy to the global velocity, see getWritableVelocity for a detailed explanation
 */
inline void applyWriteableVelocity(const std::vector<distributed_vector>& writeable, std::vector<distributed_vector>& member){
	TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
	assert (member.size() == writeable.size());
	for (size_t i = 0; i < member.size(); i++){
		member = writeable;
	}
}

/**
 * @short get a writeable copy of the velocity
 */
inline distributed_vector& getWriteableDensity(distributed_vector& writeable, const distributed_vector& member, const dealii::IndexSet& locally_owned){
	TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
	writeable.reinit(locally_owned, member.get_mpi_communicator(), true);
	assert (not writeable.has_ghost_elements());
	writeable = member;
	assert (not writeable.has_ghost_elements());
	return writeable;
}

/**
 * @short copy the changes in the writeable copy to the global density, see getWritableVelocity for a detailed explanation
 */
inline void applyWriteableDensity(const distributed_vector& writeable, distributed_vector& member){
	TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
	member = writeable;
}


} /* CFDSolverUtilities */
} /* namespace natrium */

#endif /* CFDSOLVERUTILITIES_H_ */
