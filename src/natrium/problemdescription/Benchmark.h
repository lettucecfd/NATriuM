/**
 * @file Benchmark.h
 * @short Abstract class for the description of a CFD Benchmark problem.
 * @date 26.05.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */


#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include "problemdescription/ProblemDescription.h"

#include "deal.II/grid/tria.h"

#include "BoundaryCollection.h"
#include "../utilities/BasicNames.h"


namespace natrium {

template<size_t dim> class Benchmark: public ProblemDescription<dim> {
public:
	Benchmark(
			shared_ptr<dealii::Triangulation<dim> > triangulation,
			double viscosity, double characteristicLength) ;
	virtual ~Benchmark(){
	}

	/**
	 * @short get Analytic velocity at one point in space and time
	 */
	virtual void getAnalyticVelocity(const dealii::Point<dim>& x, double t, dealii::Point<dim>& velocity) const = 0;


	/**
	 * @short get full analytic solution
	 */
	void getAnalyticSolution(double time, vector<distributed_vector>& analyticSolution,
			const vector<dealii::Point<2> >& supportPoints) const {
		assert(analyticSolution.at(0).size() == supportPoints.size());
		assert(analyticSolution.at(1).size() == supportPoints.size());
		if (dim ==3)
			assert(analyticSolution.at(2).size() == supportPoints.size());
		assert(supportPoints.size() > 0);
		dealii::Point<dim> velocity;
		for (size_t i = 0; i < supportPoints.size(); i++) {
			getAnalyticVelocity(supportPoints.at(i),
								time, velocity);
			analyticSolution.at(0)(i) = velocity(0);
			analyticSolution.at(1)(i) = velocity(1);
			/*if (dim == 3){
				analyticSolution.at(2)(i) == velocity(2);
			}*/
		}
	}
};

template<size_t dim>
inline Benchmark<dim>::Benchmark(
		shared_ptr<dealii::Triangulation<dim> > triangulation,
		double viscosity, double characteristicLength): ProblemDescription<dim>(
				triangulation, viscosity, characteristicLength)
				{
}

} /* namespace natrium */

#endif /* BENCHMARK_H_ */
