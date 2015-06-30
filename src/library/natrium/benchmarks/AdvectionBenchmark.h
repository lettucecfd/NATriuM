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

void getAnalyticSolution(double time, distributed_vector& analyticSolution,
		const vector<dealii::Point<2> >& supportPoints, bool is_smooth = true);



AdvectionResult oneTest(size_t refinementLevel, size_t fe_order, double deltaT,
		size_t numberOfTimeSteps, bool is_smooth = true,  bool output_to_std_dir = false, bool useCentralFlux = false);


} /*namespace  AdvectionBenchmark */
} /* namespace natrium */

#endif /* ADVECTIONBENCHMARK_H_ */
