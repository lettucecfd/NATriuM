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
protected:
	/// function to define initial densities
	shared_ptr<dealii::Function<dim> > m_analyticRho;

	/// function to define initial velocities (u,v,w for x,y,z component)
	shared_ptr<dealii::Function<dim> > m_analyticU;
	shared_ptr<dealii::Function<dim> > m_analyticV;
	shared_ptr<dealii::Function<dim> > m_analyticW;
public:
	/// Constructor
	Benchmark(shared_ptr<Mesh<dim> > triangulation,
			double viscosity, double characteristicLength);

	/// Destructor
	virtual ~Benchmark() {
	}

	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getAnalyticRhoFunction(double time) const {
		m_analyticRho->set_time(time);
		return m_analyticRho;
	}

	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getAnalyticUFunction(double time) const {
		m_analyticU->set_time(time);
		return m_analyticU;
	}


	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getAnalyticVFunction(double time) const {
		m_analyticV->set_time(time);
		return m_analyticV;
	}


	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getAnalyticWFunction(double time) const {
		m_analyticW->set_time(time);
		return m_analyticW;
	}


	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getInitialRhoFunction() const {
		m_analyticRho->set_time(0);
		return m_analyticRho;
	}


	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getInitialUFunction() const {
		m_analyticU->set_time(0);
		return m_analyticU;
	}


	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getInitialVFunction() const {
		m_analyticV->set_time(0);
		return m_analyticV;
	}


	/*
	 * @short analytic solution at time t
	 * @note This function sets the time at which the analytic solution is evaluated.
	 * It is not recommended to store the return value in any local object that is used
	 * in more than one time step.
	 */
	const shared_ptr<dealii::Function<dim> >& getInitialWFunction() const {
		m_analyticW->set_time(0);
		return m_analyticW;
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
