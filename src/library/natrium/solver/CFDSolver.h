/**
 * @file CFDSolver.h
 * @short Central class of the CFD Simulation based on the Discrete Boltzmann Equation (DBE)
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVER_H_
#define CFDSOLVER_H_

#include <exception>

#include "deal.II/numerics/data_out.h"

#include "SolverConfiguration.h"
#include "DistributionFunctions.h"
#include "SolverStats.h"
#include "TurbulenceStats.h"

#include "../dataprocessors/DataProcessor.h"

#include "../problemdescription/ProblemDescription.h"

#include "../advection/AdvectionOperator.h"
#include "../advection/SEDGMinLee.h"

#include "../collision/CollisionModel.h"

#include "../smoothing/Filter.h"

#include "../timeintegration/TimeIntegrator.h"

#include "../utilities/BasicNames.h"
#include "../utilities/Math.h"
#include "../utilities/NATriuMException.h"
#include "../utilities/Timing.h"

namespace natrium {

/* forward declarations */
class Stencil;

/**
 * @short Exception class for CFDSolver
 */
class CFDSolverException: public NATriuMException {
private:
	std::string message;
public:
	CFDSolverException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	CFDSolverException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~CFDSolverException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};


/** @short The central class for the CFD simulation based on the DBE.
 *  @note  The CFDSolver itself is quite static but it contains interchangeable modules, e.g. for the
 *         Stencil or the time integrator. By these means, a variety of different simulation
 *         methods can be covered.
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class CFDSolver {
	template<size_t dim2> friend class SolverStats;
	template<size_t dim3> friend class TurbulenceStats;

private:
	/// particle distribution functions
	DistributionFunctions m_f;

	/// macroscopic density
	distributed_vector m_density;

	/// macroscopic velocity
	vector<distributed_vector> m_velocity;

	/// macroscopic density
	distributed_vector m_tmpDensity;

	/// temporary velocity to test for convergence
	vector<distributed_vector> m_tmpVelocity;

	/// description of the CFD problem (boundraries, initial values, etc.)
	boost::shared_ptr<ProblemDescription<dim> > m_problemDescription;

	/// global streaming data
	boost::shared_ptr<AdvectionOperator<dim> > m_advectionOperator;

	/// DdQq Boltzmann model (e.g. D2Q9)
	boost::shared_ptr<Stencil> m_stencil;

	/// Description of the collision algorithm
	boost::shared_ptr<CollisionModel> m_collisionModel;

	/// Time Integrator for the solution of the ODE, which stems from the space discretization
	boost::shared_ptr<
			TimeIntegrator<distributed_sparse_block_matrix,
					distributed_block_vector> > m_timeIntegrator;

	/// Configuration of the solver
	boost::shared_ptr<SolverConfiguration> m_configuration;

	/// Filtering scheme
	boost::shared_ptr<Filter<dim> > m_filter;

	/// the number of the first iteration (normally 0, except for restart at a checkpoint)
	size_t m_iterationStart;

	/// the physical time passed (normally initialized with 0.0, except for restart at a checkpoint)
	double m_time;

	/// a vector that indicates if a dofs is at the boundary (for each dof)
	vector<bool> m_isDoFAtBoundary;

	/// current iteration
	size_t m_i;

	/// table out
	boost::shared_ptr<SolverStats<dim> > m_solverStats;
	boost::shared_ptr<TurbulenceStats<dim> > m_turbulenceStats;

	// vector of data processors
	vector< boost::shared_ptr<DataProcessor<dim> > > m_dataProcessors;

	// starting time
	time_t m_tstart;

	// residuum
	double m_residuumDensity;
	double m_residuumVelocity;

	// vector for nonlinear boundary conditions
	distributed_block_vector m_boundaryVector;

	// vector of grid points
	map<dealii::types::global_dof_index, dealii::Point<dim> > m_supportPoints;




protected:

	/// save the distribution functions to files for checkpointing
	void saveDistributionFunctionsToFiles(const string& directory);

	/// load the distribution functions from files for checkpointing
	void loadDistributionFunctionsFromFiles(const string& directory);

	/// gives the possibility for Benchmark instances to add the analytic solution to output
	virtual void addAnalyticSolutionToOutput(dealii::DataOut<dim>&) {
	}

public:

	/// constructor
	/// @note: has to be inlined, if the template parameter is not made explicit
	CFDSolver(boost::shared_ptr<SolverConfiguration> configuration,
			boost::shared_ptr<ProblemDescription<dim> > problemDescription);

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
	 * @short filter solution
	 */
	void filter();

	/**
	 * @short test for stop conditions
	 */
	bool stopConditionMet();

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	void applyInitialDensities(distributed_vector& initialDensities,
			const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const;

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	void applyInitialVelocities(vector<distributed_vector>& initialVelocities,
			const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const;

	/**
	 * @short create output data and write to file
	 */
	virtual void output(size_t iteration);

	/**
	 *
	 */
	bool hasGeometryChanged() {
		return false;
	}

	const distributed_vector& getDensity() const {
		return m_density;
	}

	const DistributionFunctions& getF() const {
		return m_f;
	}

	const vector<distributed_vector>& getVelocity() const {
		return m_velocity;
	}

	const boost::shared_ptr<AdvectionOperator<dim> >& getAdvectionOperator() const {
		return m_advectionOperator;
	}

	const boost::shared_ptr<Stencil>& getStencil() const {
		return m_stencil;
	}

	const boost::shared_ptr<CollisionModel>& getCollisionModel() const {
		return m_collisionModel;
	}

	const boost::shared_ptr<SolverConfiguration>& getConfiguration() const {
		return m_configuration;
	}

	const boost::shared_ptr<ProblemDescription<dim> >& getProblemDescription() const {
		return m_problemDescription;
	}

	const boost::shared_ptr<
			TimeIntegrator<distributed_vector, distributed_sparse_matrix> >& getTimeIntegrator() const {
		return m_timeIntegrator;
	}

	size_t getNumberOfDoFs() const {
		return m_advectionOperator->getNumberOfDoFs();
	}
	double getMaxVelocityNorm() const {
		double max = m_velocity.at(0).linfty_norm();
		double comp2 = m_velocity.at(1).linfty_norm();
		if (comp2 > max){
			max = comp2;
		}
		double comp3 = 0;
		if (dim == 3)
			comp3 = m_velocity.at(2).linfty_norm();
		if (comp3 > max){
			max = comp3;
		}
		return max;
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

	double getTime() const {
		return m_time;
	}

	size_t getIteration() const {
		return m_i;
	}

	void setIteration(size_t iteration) {
		m_i = iteration;
	}

	const boost::shared_ptr<SolverStats<dim> >& getSolverStats() const {
		return m_solverStats;
	}

	double getTau() const;

	double getResiduumDensity() const {
		return m_residuumDensity;
	}

	double getResiduumVelocity() const {
		return m_residuumVelocity;
	}

	void printRuntimeSummary() const {
		pout << Timing::getOutStream().str() << endl;
	}

	void addToVelocity(boost::shared_ptr<dealii::Function<dim> > function);

	void scaleVelocity(double scaling_factor);

	const map<dealii::types::global_dof_index, dealii::Point<dim> >& getSupportPoints() const {
		return m_supportPoints;
	}

	void setVelocity(const vector<distributed_vector>& velocity) {
		m_velocity = velocity;
	}

	void appendDataProcessor( boost::shared_ptr<DataProcessor<dim> > proc){
		m_dataProcessors.push_back(proc);
	}
}
;

} /* namespace natrium */

#endif /* CFDSOLVER_H_ */
