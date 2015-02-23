/*
 * DealIIWrapper_test.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: kraemer
 */

#include <timeintegration/DealIIWrapper.h>

#include "boost/test/unit_test.hpp"

#include "deal.II/lac/sparsity_pattern.h"
#include "deal.II/lac/compressed_sparsity_pattern.h"

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(DealIIWrapper_test)

BOOST_AUTO_TEST_CASE(DealIIWrapper_Convergence_test) {
	cout << "DealIIWrapper_Convergence_test..." << endl;

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
		DealIIWrapper<distributed_sparse_matrix, distributed_vector> tm(dt);
		distributed_vector f(1);
		f(0) = 1;
		double t = 0;
		for (size_t i = 0; i < numberOfSteps; i++) {
			t += dt;
			tm.step(f, A, g);
			double error = fabs(f(0) - exp(lambda * t));
			BOOST_ASSERT(error < 1e-3);
		}
	}
	cout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(DealIIWrapper_MultiBlock_test) {
	cout << "DealIIWrapper_MultiBlock_test..." << endl;

	// solve ODE F(f) = A*f
	// with initial value f0 = [ 1 2 ]
	// and A = [ 1 -1 , 0 3 ]
	// analytic solution: f(x) = e^(A*x)
	const double dt = 0.001;
	const size_t numberOfSteps = 1000;
	// build matrix
	distributed_sparse_block_matrix A;
#ifdef WITH_TRILINOS
	dealii::CompressedSparsityPattern compressedSparsityPattern(1, 1);
	compressedSparsityPattern.add(0, 0);
	A.reinit(2, 2);
	A.block(0, 0).reinit(compressedSparsityPattern);
	A.block(0, 1).reinit(A.block(0, 0));
	A.block(1, 0).reinit(A.block(0, 0));
	A.block(1, 1).reinit(A.block(0, 0));
	A.collect_sizes();
#else
	dealii::BlockCompressedSparsityPattern cSparse(2, 2);
	for (size_t I = 0; I < 2; I++) {
		for (size_t J = 0; J < 2; J++) {
			cSparse.block(I, J).reinit(1, 1);
		}
	}
	cSparse.collect_sizes();
	dealii::BlockSparsityPattern sparse(2, 2);
	for (size_t I = 0; I < 2; I++) {
		for (size_t J = 0; J < 2; J++) {
			sparse.block(I, J).copy_from(cSparse.block(I, J));
		}
	}
	sparse.collect_sizes();
	A.reinit(sparse);
#endif
	A.set(0, 0, 1);
	A.set(0, 1, -1);
	A.set(1, 0, 0);
	A.set(1, 1, 3);

	for (double theta = 0.0; theta <= 1; theta += 0.5) {
		DealIIWrapper<distributed_sparse_block_matrix, distributed_block_vector> tm(
				dt);
		// initialize block vectors
		distributed_block_vector f;
		distributed_block_vector b;
#ifdef WITH_TRILINOS
		f.reinit(2);
		for (size_t i = 0; i < 2; i++) {
			f.block(i).reinit(1);
		}
		f.collect_sizes();
		//b
		b.reinit(2);
		for (size_t i = 0; i < 2; i++) {
			b.block(i).reinit(1);
		}
		b.collect_sizes();
#else
		f.reinit(2, 1);
		b.reinit(2, 1);
#endif
		f(0) = 1;
		f(1) = 2;
		b(0) = 0;
		b(1) = 0;
		double c0 = f(0);
		double c1 = f(1);
		double t = 0;
		for (size_t i = 0; i < numberOfSteps; i++) {
			t += dt;
			tm.step(f, A, b);
		}
		// calculate analytic solution
		// f1(x) = c1 * exp(x) - 0.5 * c2 * exp(x) [ exp(2x) - 1]
		// f2(x) = c2 * exp(3x)
		// with c2 = f2(0) and c1 = f1(0)
		double f0 = c0 * exp(t) - 0.5 * c1 * exp(t) * (exp(2*t) - 1.0);
		double f1 = c1 * exp(3*t);
		// compare
		BOOST_CHECK_CLOSE( (double) f(0), f0, 1e-1);
		BOOST_CHECK_CLOSE( (double) f(1), f1, 1e-1);
	}
	A.clear();
	cout << "done." << endl;
} /* DealIIWrapper_MultiBlock_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
