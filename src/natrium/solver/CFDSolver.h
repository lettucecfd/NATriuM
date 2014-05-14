/**
 * @file CFDSolver.h
 * @short Central class of the CFD Simulation based on the Discrete Boltzmann Equation (DBE)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVER_H_
#define CFDSOLVER_H_

#include <exception>

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
#include "DistributionFunctions.h"

namespace natrium {

/**
 * @short Exception class for CFDSolver
 */
class CFDSolverException: public std::exception {
private:
	std::string message;
public:
	CFDSolverException(const char *msg) :
			message(msg) {
	}
	~CFDSolverException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/** @short The central class for the CFD simulation based on the DBE.
 *  @note  The CFDSolver itself is quite static but it contains interchangeable modules, e.g. for the
 *         Boltzmann model or the time integrator. By these means, a variety of different simulation
 *         methods can be covered.
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class CFDSolver {

private:

	/// particle distribution functions
	DistributionFunctions m_f;

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
	shared_ptr<TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector> > m_timeIntegrator;

	/// Configuration of the solver
	shared_ptr<SolverConfiguration> m_configuration;

	/// the number of the first iteration (normally 0, except for restart at a checkpoint)
	size_t m_iterationStart;

	/// a vector that indicates if a dofs is at the boundary (for each dof)
	vector<bool> m_isBoundary;

protected:

	/// save the distribution functions to files for checkpointing
	void saveDistributionFunctionsToFiles(const string& directory);

	/// load the distribution functions from files for checkpointing
	void loadDistributionFunctionsFromFiles(const string& directory);



public:

	/// constructor
	/// @note: has to be inlined, if the template parameter is not made explicit
	CFDSolver(shared_ptr<SolverConfiguration> configuration,
			shared_ptr<ProblemDescription<dim> > problemDescription);

/// destructor
	virtual ~CFDSolver() {
	}
	;

	/**
	 * @short initialize distribution functions
	 */
	void initializeDistributions();

	/**
	 * @short Advection in all directions
	 */
	void stream();

	/**
	 *  @short Low-level collide function
	 */
	void collide();

	/**
	 * @short reassembly of all matrices
	 */
	void reassemble();

	/**
	 * @short run CFD solver
	 */
	void run();

	/**
	 * @short create output data and write to file
	 */
	void output(size_t iteration);

	const distributed_vector& getDensity() const {
		return m_density;
	}

	const vector<distributed_vector>& getF() const {
		return m_f;
	}

	const vector<distributed_vector>& getVelocity() const {
		return m_velocity;
	}

	const shared_ptr<AdvectionOperator<dim> >& getAdvectionOperator() const {
		return m_advectionOperator;
	}

	const shared_ptr<BoltzmannModel>& getBoltzmannModel() const {
		return m_boltzmannModel;
	}

	const shared_ptr<CollisionModel>& getCollisionModel() const {
		return m_collisionModel;
	}

	const shared_ptr<SolverConfiguration>& getConfiguration() const {
		return m_configuration;
	}

	const shared_ptr<ProblemDescription<dim> >& getProblemDescription() const {
		return m_problemDescription;
	}

	const shared_ptr<TimeIntegrator<distributed_vector, distributed_sparse_matrix> >& getTimeIntegrator() const {
		return m_timeIntegrator;
	}

	size_t getNumberOfDoFs() const {
		return m_advectionOperator->getNumberOfDoFs();
	}

	double getMaxVelocityNorm() const {
		double maxnorm = 0.0;
		for (size_t i = 0; i < getNumberOfDoFs(); i++) {
			double norm = 0.0;
			for (size_t j = 0; j < dim; j++) {
				norm += m_velocity.at(j)(i) * m_velocity.at(j)(i);
			}
			if (norm > maxnorm) {
				maxnorm = norm;
			}
		}
		return sqrt(maxnorm);
	}

	double getMaxDensityDeviationFrom(double referenceDensity) const {
		double maxdev = 0.0;
		for (size_t i = 0; i < getNumberOfDoFs(); i++) {
			double dev = fabs(m_density(i) - referenceDensity);
			if (dev > maxdev) {
				maxdev = dev;
			}
		}
		return maxdev;
	}

	size_t getIterationStart() const {
		return m_iterationStart;
	}
}
;

} /* namespace natrium */


#endif /* CFDSOLVER_H_ */
