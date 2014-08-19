/**
 * @file ThetaMethod_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "timeintegration/ThetaMethod.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/lac/sparsity_pattern.h"
#include "deal.II/lac/compressed_sparsity_pattern.h"

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(ThetaMethod_test)

BOOST_AUTO_TEST_CASE(ThetaMethod_Convergence_test) {
	cout << "ThetaMethod_Convergence_test..." << endl;

	// solve ODE F(f) = lambda*f
	// with initial value f0 = 1
	// analytic solution: f(x) = e^(lamda*x)
	const double dt = 0.001;
	const size_t numberOfSteps = 1000;
	const double lambda = -2.0;
	distributed_vector g(1);
	g(0) = 0;
	// build the 1x1 matrix [[lambda]]
	dealii::CompressedSparsityPattern compressedSparsityPattern(1, 1);
	compressedSparsityPattern.add(0, 0);
	dealii::SparsityPattern sparsityPattern;
	sparsityPattern.copy_from(compressedSparsityPattern);
	distributed_sparse_matrix A;
	A.reinit(sparsityPattern);
	A.set(0, 0, lambda);

	for (double theta = 0.0; theta <= 1; theta += 0.5) {
		ThetaMethod<distributed_sparse_matrix, distributed_vector> tm(dt, 1,
				theta);
		distributed_vector f(1);
		f(0) = 1;
		double t = 0;
		for (size_t i = 0; i < numberOfSteps; i++) {
			t += dt;
			tm.step(f, A, g);
			double error = fabs(f(0) - exp(lambda * t));
			cout << f(0) << " " << exp(lambda * t) << endl;
			BOOST_ASSERT(error < 1e-3);
		}
	}
	cout << "done." << endl;
}
BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
