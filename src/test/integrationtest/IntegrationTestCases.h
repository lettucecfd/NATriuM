/*
 * IntegrationTestCases.h
 *
 *  Created on: Jan 28, 2015
 *      Author: kraemer
 */

#ifndef INTEGRATIONTESTCASES_H_
#define INTEGRATIONTESTCASES_H_

#include <sstream>

#include "natrium/utilities/BasicNames.h"

namespace natrium {
namespace IntegrationTestCases {

struct TestResult {
	size_t id;
	string name;
	string details;
	vector<string> quantity;
	vector<double> expected;
	vector<double> threshold;
	vector<double> outcome;
	shared_ptr<std::stringstream> error_msg;
	double time;
	bool success;
	TestResult(){
		id = 0;
		time = 0;
		success = 0;
		error_msg = make_shared<std::stringstream>();
	};
	~TestResult(){};
};

TestResult ConvergenceTestPeriodic ();

TestResult ConvergenceTestMovingWall ();

TestResult ConvergenceTestImplicitLBM ();

TestResult ConvergenceTestExponentialLBM();


} /* namespace IntegrationTests */
} /* namespace natrium */

#endif /* INTEGRATIONTESTCASE_H_ */
