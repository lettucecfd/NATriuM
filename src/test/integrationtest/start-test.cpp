/**
 * @file start-test.cpp
 * @short Execute integration test cases
 * @date 29.01.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "IntegrationTestCases.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/HtmlTrace.h"
#include "natrium/utilities/MPIGuard.h"

using namespace natrium;

void print_line_html(IntegrationTestCases::TestResult& result,
		std::ofstream& html) {

	if (is_MPI_rank_0()) {
//////////////////////
// Print html trace //
//////////////////////
		size_t n = result.quantity.size();
		assert(n > 0);

		//begin line
		html << "  <tr align=\"center\">\n";
		html << "    <td rowspan='" << n << "'>" << result.id << "</td>\n";
		html << "    <td rowspan='" << n << "'>" << result.name << "</td>\n";

		// quantity name (first measured qty)
		html << "    <td>" << result.quantity.at(0) << "</td>\n";
		// expected value
		html << "    <td>" << result.expected.at(0) << "</td>\n";
		// simulated value
		if (fabs(result.outcome.at(0) - result.expected.at(0))
				<= result.threshold.at(0)) {
			html << "    <td bgcolor='LightGreen'>" << result.outcome.at(0)
					<< "</td>\n";
		} else {
			html << "    <td bgcolor='red'>" << result.outcome.at(0)
					<< "</td>\n";
		}
		// success
		if (result.success) {
			html << "    <td rowspan='" << n
					<< "' bgcolor='LightGreen' align='center'> + </td>\n";
		} else {
			html << "    <td rowspan='" << n
					<< "' bgcolor='red' align='center'> - </td>\n";
		}
		// time
		html << "    <td rowspan='" << n << "'> " << result.time << "</td>\n";
		// details
		html << "    <td rowspan='" << n << "'> " << result.details
				<< "</td>\n";
		// error details
		if (result.success) {
			html << "    <td rowspan='" << n << "'> </td>\n";
		} else {
			html << "    <td rowspan='" << n << "' bgcolor='red'> "
					<< result.error_msg->str().c_str() << "</td>\n";
		}
		//end line
		html << "  </tr>\n";

		// other measured quantities
		for (size_t i = 1; i < n; i++) {
			//begin line
			html << "  <tr align=\"center\">\n";
			// quantity name
			html << "    <td>" << result.quantity.at(i) << "</td>\n";
			// expected value
			html << "    <td>" << result.expected.at(i) << "</td>\n";
			// simulated value
			if (fabs(result.outcome.at(i) - result.expected.at(i))
					<= result.threshold.at(i)) {
				html << "    <td bgcolor='LightGreen'>" << result.outcome.at(i)
						<< "</td>\n";
			} else {
				html << "    <td bgcolor='red'>" << result.outcome.at(i)
						<< "</td>\n";
			}
			//end line
			html << "  </tr>\n";
		}
		html.flush();
	} /*	if (0 == dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) */

}

bool print_result(IntegrationTestCases::TestResult& result, std::ofstream & out, bool errors){
	bool err = errors;
	print_line_html(result, out);
	if (result.success) {
		pout << "-  " << result.name << " ... " << "OK." << endl;
	} else {
		pout << "-  " << result.name << " ... " << "Error: "
				<< result.error_msg->str().c_str()
				<< " See natrium.html for details." << endl;
		err = true;
	}
	return err;
}


// Main function
int main(int argc, char **argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Start integration tests. This will take a few minutes..." << endl;

	bool errors = false;
	HtmlTrace htmlTrace(is_MPI_rank_0());
	IntegrationTestCases::TestResult result;

	// Test 15: Convergence Moving Walls (SL)
	result = IntegrationTestCases::ConvergenceTestMovingWall(SEMI_LAGRANGIAN);
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 14: Convergence of all collision schemes
	result = IntegrationTestCases::ConvergenceTestCollisionSchemes();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 13: Convergence of semi-Lagrangian advection scheme (nonsmooth problem)
	result = IntegrationTestCases::ConvergenceTestSemiLagrangianAdvectionNonsmooth();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 12: Convergence of semi-Lagrangian advection scheme (smooth problem)
	result = IntegrationTestCases::ConvergenceTestSemiLagrangianAdvectionSmooth();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 11: Convergence of semi-Lagrangian for periodic flow
	result = IntegrationTestCases::ConvergenceTestSemiLagrangianPeriodic();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 10: Convergence Forcing Schemes 3D
	result = IntegrationTestCases::ConvergenceTestForcingSchemes3D();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 9: Convergence Forcing Schemes 2D
	result = IntegrationTestCases::ConvergenceTestForcingSchemes2D();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 8: Convergence Moving Walls
	result = IntegrationTestCases::ConvergenceTestMovingWall(SEDG);
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 7: Convergence 3D flow
	result = IntegrationTestCases::ConvergenceTest3D();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 6: Convergence DealIIWrapper
	result = IntegrationTestCases::ConvergenceTestDealIIWrapper();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 5: Convergence Exponential LBM
	result = IntegrationTestCases::ConvergenceTestExponentialLBM();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 4: Convergence Implicit LBM
	result = IntegrationTestCases::ConvergenceTestImplicitLBM();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 3: Convergence Periodic Walls
	result = IntegrationTestCases::ConvergenceTestPeriodic();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 2: Convergence Pure Linear Advection (non-smooth)
	result = IntegrationTestCases::ConvergenceSEDGLinearAdvectionNonsmooth();
	errors = print_result(result, htmlTrace.getHtml(), errors);

	// Test 1: Convergence Pure Linear Advection (smooth)
	result = IntegrationTestCases::ConvergenceSEDGLinearAdvectionSmooth();
	errors = print_result(result, htmlTrace.getHtml(), errors);


	// FINALIZE
	if (errors) {
		pout << "Done. Errors occured in tests. See natrium.html for details."
				<< endl;
	} else {
		pout << "All tests passed. No errors detected." << endl;
	}
	return 0;
}
