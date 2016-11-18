/**
 * @file ThetaMethod_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/timeintegration/ThetaMethod.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/lac/sparsity_pattern.h"
#include "deal.II/lac/compressed_sparsity_pattern.h"

#include "natrium/benchmarks/AdvectionBenchmark.h"
#include "natrium/utilities/BasicNames.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(ThetaMethod_test)

// This unit test module ensures that the time integrator work
// for the vector classes numeric_vector, block_vector and distributed_vector.
// A test for distributed_block_vector is provided by the integration test.

BOOST_AUTO_TEST_CASE(ThetaMethod_Convergence_test) {
	pout << "ThetaMethod_Convergence_test..." << endl;

	// solve ODE F(f) = lambda*f
	// with initial value f0 = 1
	// analytic solution: f(x) = e^(lamda*x)
	const double dt = 0.001;
	const size_t numberOfSteps = 1000;
	const double lambda = -2.0;
	numeric_vector g(1);
	g(0) = 0;
	// build the 1x1 matrix [[lambda]]
	dealii::DynamicSparsityPattern compressedSparsityPattern(1, 1);
	compressedSparsityPattern.add(0, 0);
	dealii::SparsityPattern sparsityPattern;
	sparsityPattern.copy_from(compressedSparsityPattern);
	sparse_matrix A;
	A.reinit(sparsityPattern);
	A.set(0, 0, lambda);

	for (double theta = 0.0; theta <= 1; theta += 0.5) {
		ThetaMethod<sparse_matrix, numeric_vector> tm(dt, g, theta);
		numeric_vector f(1);
		f(0) = 1;
		double t = 0;
		for (size_t i = 0; i < numberOfSteps; i++) {
			t += dt;
			tm.step(f, A, g);
			double error = fabs(f(0) - exp(lambda * t));
			BOOST_ASSERT(error < 1e-3);
		}
	}
	pout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(ThetaMethod_MultiBlock_test) {
	pout << "ThetaMethod_MultiBlock_test..." << endl;

	// solve ODE F(f) = A*f
	// with initial value f0 = [ 1 2 ]
	// and A = [ 1 -1 , 0 3 ]
	// analytic solution: f(x) = e^(A*x)
	const double dt = 0.0001;
	const size_t numberOfSteps = 1000;
	// build matrix
	sparse_block_matrix A;

	dealii::BlockDynamicSparsityPattern cSparse(2, 2);
	for (size_t iI = 0; iI < 2; iI++) {
		for (size_t J = 0; J < 2; J++) {
			cSparse.block(iI, J).reinit(1, 1);
		}
	}
	cSparse.collect_sizes();
	dealii::BlockSparsityPattern sparse(2, 2);
	for (size_t iI = 0; iI < 2; iI++) {
		for (size_t J = 0; J < 2; J++) {
			sparse.block(iI, J).copy_from(cSparse.block(iI, J));
		}
	}
	sparse.collect_sizes();
	A.reinit(sparse);

	A.set(0, 0, 1);
	A.set(0, 1, -1);
	A.set(1, 0, 0);
	A.set(1, 1, 3);

	for (double theta = 0.0; theta <= 1; theta += 0.5) {
		// initialize block vectors
		block_vector f(2, 1);
		block_vector b(2, 1);
		f(0) = 1;
		f(1) = 2;
		b(0) = 0;
		b(1) = 0;
		double c0 = f(0);
		double c1 = f(1);
		double t = 0;
		ThetaMethod<sparse_block_matrix, block_vector> tm(dt, f, theta);
		for (size_t i = 0; i < numberOfSteps; i++) {
			t += dt;
			tm.step(f, A, b);
		}
		// calculate analytic solution
		// f1(x) = c1 * exp(x) - 0.5 * c2 * exp(x) [ exp(2x) - 1]
		// f2(x) = c2 * exp(3x)
		// with c2 = f2(0) and c1 = f1(0)
		double f0 = c0 * exp(t) - 0.5 * c1 * exp(t) * (exp(2 * t) - 1.0);
		double f1 = c1 * exp(3 * t);
		// compare
		BOOST_CHECK_CLOSE((double ) f(0), f0, 1e-1);
		BOOST_CHECK_CLOSE((double ) f(1), f1, 1e-1);
	}
	A.clear();
	pout << "done." << endl;
} /* ThetaMethod_MultiBlock_test */

BOOST_AUTO_TEST_CASE(ThetaMethod_MPI_test) {
	pout << "ThetaMethod_MPI_test..." << endl;

	// Advection benchmark

	size_t N = 2;	// refinement level
	size_t p = 4;	// fe order
	bool is_smooth = true;

	double delta_x = 1. / (pow(2, N));
	double delta_t;

	delta_t = 0.2 * pow(0.5, N) / ((p + 1) * (p + 1));

	double t_end = 0.1;
	AdvectionBenchmark::AdvectionResult advectionResult =
			AdvectionBenchmark::oneTest(N, p, delta_t, t_end,
					THETA_METHOD, NONE, is_smooth, false, false);

	double result = std::log10(advectionResult.normSup);
	// taken from the integration test
	double expected = std::log10(400 * std::pow(0.3 * delta_x, p + 1));

	BOOST_CHECK_SMALL(result - expected, 0.8);

	pout << "done." << endl;
}/* ThetaMethod_MPI_test */

BOOST_AUTO_TEST_SUITE_END()

