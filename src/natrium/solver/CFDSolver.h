/**
 * @file CFDSolver.h
 * @short Central class of the CFD Simulation based on the Discrete Boltzmann Equation (DBE)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVER_H_
#define CFDSOLVER_H_

#include "../problemdescription/ProblemDescription.h"

#include "../advection/AdvectionOperator.h"
#include "../advection/SEDGMinLee.h"

#include "../boltzmannmodels/BoltzmannModel.h"
#include "../boltzmannmodels/D2Q9IncompressibleModel.h"

#include "../collisionmodels/CollisionModel.h"
#include "../collisionmodels/BGKTransformed.h"

#include "../timeintegration/TimeIntegrator.h"
#include "../timeintegration/RungeKutta5LowStorage.h"

#include "../utilities/BasicNames.h"

#include "SolverConfiguration.h"

namespace natrium {

/** @short The central class for the CFD simulation based on the DBE.
 *  @note  The CFDSolver itself is quite static but it contains interchangeable modules, e.g. for the
 *         Boltzmann model or the time integrator. By these means, a variety of different simulation
 *         methods can be covered.
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class CFDSolver {

private:

	/// particle distribution functions
	vector<distributed_vector> m_f;

	/// macroscopic density
	distributed_vector m_density;

	/// macroscopic velocity
	vector<distributed_vector> m_velocity;

	/// description of the CFD problem (boundraries, initial values, etc.)
	shared_ptr<ProblemDescription<dim> > m_problemDescription;

	/// global streaming data
	shared_ptr<AdvectionOperator<dim> > m_advectionOperator;

	/// DdQq Boltzmann model (e.g. D2Q9)
	shared_ptr<BoltzmannModel> m_boltzmannModel;

	/// Description of the collision algorithm
	shared_ptr<CollisionModel> m_collisionModel;

	/// Time Integrator for the solution of the ODE, which stems from the space discretization
	shared_ptr<TimeIntegrator> m_timeIntegrator;

	/// Configuration of the solver
	shared_ptr<SolverConfiguration> m_configuration;

public:

	/// constructor
	/// @note: has to be inlined, if the template parameter is not made explicit
	CFDSolver(shared_ptr<SolverConfiguration> configuration,
			shared_ptr<ProblemDescription<dim> > problemDescription) {

		/// check if problem and solver configuration fit together
		configuration->checkProblem(problemDescription);
		m_problemDescription = problemDescription;
		m_configuration = configuration;

		/// Build streaming data object
		if (Streaming_MinLee2011 == configuration->getStreamingDataType()) {
/*			m_streamingData = make_shared<DataMinLee2011<dim> >(
					m_problemDescription->getTriangulation(),
					configuration->getOrderOfFiniteElement());
*/		}

		/// Build boltzmann model
		if (Stencil_D2Q9 == configuration->getStencilType()) {
			m_boltzmannModel = make_shared<D2Q9IncompressibleModel>();
		}

		/// Build collision model
		if (Collision_BGKTransformed == configuration->getCollisionType()){
			m_collisionModel = make_shared<BGKTransformed>(m_problemDescription->getRelaxationParameter(), m_boltzmannModel);
		}

		/// Build time integrator
		/*if (Integrator_RungeKutta5LowStorage == configuration->getTimeIntegratorType()){
			m_timeIntegrator = make_shared<RungeKutta5LowStorage>();
		}*/

	}
	;

/// destructor
	virtual ~CFDSolver() {
	}
	;

	/**
	 * @short run CFD solver
	 */
	void run(){};

};

} /* namespace natrium */



#endif /* CFDSOLVER_H_ */
