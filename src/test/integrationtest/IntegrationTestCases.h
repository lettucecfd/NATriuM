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
	boost::shared_ptr<std::stringstream> error_msg;
	double time;
	bool success;
	TestResult(){
		id = 0;
		time = 0;
		success = 0;
		error_msg = boost::make_shared<std::stringstream>();
	};
	~TestResult(){};
};

TestResult ConvergenceSEDGLinearAdvectionSmooth();

TestResult ConvergenceSEDGLinearAdvectionNonsmooth();

TestResult ConvergenceTestPeriodic ();

TestResult ConvergenceTestImplicitLBM ();

TestResult ConvergenceTestExponentialLBM();

TestResult ConvergenceTestDealIIWrapper();

TestResult ConvergenceTest3D();

TestResult ConvergenceTestMovingWall ();

TestResult ConvergenceTestForcingSchemes2D ();

TestResult ConvergenceTestForcingSchemes3D ();

TestResult ConvergenceTestSemiLagrangianPeriodic ();

TestResult ConvergenceTestSemiLagrangianAdvectionSmooth ();

TestResult ConvergenceTestSemiLagrangianAdvectionNonsmooth ();


} /* namespace IntegrationTests */
} /* namespace natrium */

#endif /* INTEGRATIONTESTCASE_H_ */
