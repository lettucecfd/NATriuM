/*
 * AdvectionBenchmark.h
 *
 *  Created on: May 11, 2015
 *      Author: kraemer
 */

#ifndef ADVECTIONBENCHMARK_H_
#define ADVECTIONBENCHMARK_H_


#include "../utilities/BasicNames.h"

namespace natrium {

namespace AdvectionBenchmark {

struct AdvectionResult {
	size_t refinementLevel;
	size_t fe_order;
	double deltaX;
	double deltaT;
	double norm2;
	double normSup;
	double timesec;
};

/* function to calculate analytic solution at a point x*/
double analytic_solution(double time, const dealii::Point<2>& x,
		bool is_smooth);

void getAnalyticSolution(double time, distributed_vector& analyticSolution,
		const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints,
		const shared_ptr<AdvectionOperator<2> > & streaming,
		bool is_smooth = true);



AdvectionResult oneTest(size_t refinementLevel, size_t fe_order, double deltaT,
		size_t numberOfTimeSteps, bool is_smooth = true,  bool output_to_std_dir = false, bool useCentralFlux = false);


} /*namespace  AdvectionBenchmark */
} /* namespace natrium */

#endif /* ADVECTIONBENCHMARK_H_ */
