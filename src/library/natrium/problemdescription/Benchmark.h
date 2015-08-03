/**
 * @file Benchmark.h
 * @short Abstract class for the description of a CFD Benchmark problem.
 * @date 26.05.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BENCHMARK_H_
#define BENCHMARK_H_



#include "deal.II/grid/tria.h"

#include "BoundaryCollection.h"

#include "../problemdescription/ProblemDescription.h"

#include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim> class Benchmark: public ProblemDescription<dim> {
public:
	/// Constructor
	Benchmark(shared_ptr<Mesh<dim> > triangulation,
			double viscosity, double characteristicLength);

	/// Destructor
	virtual ~Benchmark() {
	}

	/**
	 * @short get Analytic velocity at one point in space and time
	 */
	virtual void getAnalyticVelocity(const dealii::Point<dim>& x, double t,
			dealii::Point<dim>& velocity) const = 0;

	/**
	 * @short get Analytic density at one point in space and time
	 */
	virtual double getAnalyticDensity(const dealii::Point<dim>& x,
			double t) const {
		return 1.0;
	}

	/**
	 * @short get full analytic solution for the velocity field at time t
	 */
	void getAllAnalyticVelocities(double time,
			vector<distributed_vector>& analyticSolution,
			const vector<dealii::Point<dim> >& supportPoints) const {
		// check dimensions
		assert(analyticSolution.at(0).size() == supportPoints.size());
		assert(analyticSolution.at(1).size() == supportPoints.size());
		if (dim == 3)
			assert(analyticSolution.at(2).size() == supportPoints.size());
		assert(supportPoints.size() > 0);

		// get analytic velocities
		dealii::Point<dim> velocity;
		for (size_t i = 0; i < supportPoints.size(); i++) {
			getAnalyticVelocity(supportPoints.at(i), time, velocity);
			analyticSolution.at(0)(i) = velocity(0);
			analyticSolution.at(1)(i) = velocity(1);
			if (dim == 3) {
				analyticSolution.at(2)(i) = velocity(2);
			}
		}
	}

	virtual void applyInitialDensities(distributed_vector& initialDensities,
			const vector<dealii::Point<dim> >& supportPoints) const {
		getAllAnalyticDensities(0.0, initialDensities, supportPoints);
	}

	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<dim> >& supportPoints) const {
		assert(
				initialVelocities.at(0).size()
						== initialVelocities.at(1).size());
		assert(initialVelocities.size() == dim);
		getAllAnalyticVelocities(0.0, initialVelocities, supportPoints);
	}

	/**
	 * @short get full analytic solution for the density field at time t
	 */
	virtual void getAllAnalyticDensities(double time,
			distributed_vector& analyticSolution,
			const vector<dealii::Point<dim> >& supportPoints) const {
		assert(analyticSolution.size() == supportPoints.size());
		for (size_t i = 0; i < analyticSolution.size(); i++) {
			analyticSolution(i) = getAnalyticDensity(supportPoints.at(i), time);
		}
	}

};

template<size_t dim>
inline Benchmark<dim>::Benchmark(
		shared_ptr<Mesh<dim> > triangulation, double viscosity,
		double characteristicLength) :
		ProblemDescription<dim>(triangulation, viscosity, characteristicLength) {
}

} /* namespace natrium */

#endif /* BENCHMARK_H_ */
