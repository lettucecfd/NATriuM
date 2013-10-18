/**
 * @file CFDSolver.h
 * @short Central class of the CFD Simulation based on the Discrete Boltzmann Equation (DBE)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVER_H_
#define CFDSOLVER_H_


#include "boost/shared_ptr.hpp"

#include "../streamingdata/StreamingData.h"
#include "../problemdescription/ProblemDescription.h"
#include "../boltzmannmodels/BoltzmannModel.h"
#include "../collisionmodels/CollisionModel.h"
#include "../timeintegration/TimeIntegrator.h"

using boost::shared_ptr;

namespace natrium {

/** @short The central class for the CFD simulation based on the DBE.
 *  @note  The CFDSolver itself is quite static but it contains interchangeable modules, e.g. for the
 *         Boltzmann model or the time integrator. By these means, a variety of different simulation
 *         methods can be covered.
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<int dim> class CFDSolver {

private:

	/// global streaming data
	shared_ptr<StreamingData<dim> > streamingData;

	/// description of the CFD problem (boundraries, initial values, etc.)
	shared_ptr<ProblemDescription<dim> > problemDescription;

	/// DdQq Boltzmann model (e.g. D2Q9)
	shared_ptr<BoltzmannModel> boltzmannModel;

	/// Description of the collision algorithm
	shared_ptr<CollisionModel> collisionModel;

	/// Time Integrator for the solution of the ODE, which stems from the space discretization
	shared_ptr<TimeIntegrator> timeIntegrator;


public:

	/// constructor
	/// @note: has to be inlined, if the template parameter is not made explicit
	CFDSolver(){};

	/// destructor
	virtual ~CFDSolver(){};

};

} /* namespace natrium */
#endif /* CFDSOLVER_H_ */
