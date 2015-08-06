/*
 * DealIIWrapper_test.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: kraemer
 */

#include "natrium/timeintegration/DealIIWrapper.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/lac/sparsity_pattern.h"
#include "deal.II/lac/compressed_sparsity_pattern.h"

#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(DealIIWrapper_test)

BOOST_AUTO_TEST_CASE(DealIIWrapper_Convergence_test) {
	cout << "DealIIWrapper_Convergence_test..." << endl;

	for (int n = 0; n < 12; n++) {
		DealIntegratorName name = static_cast<DealIntegratorName>(n);
		cout << "  - Integrator " << n << "...";
		// solve ODE F(f) = lambda*f
		// with initial value f0 = 1
		// analytic solution: f(x) = e^(lambda*x)
		double tmax = 1;
		double dt = 0.001;
		const double lambda = -2.0;
		UNDISTRIBUTED_VECTOR(g,1);
		g(0) = 0;
		// build the 1x1 matrix [[lambda]]
		dealii::DynamicSparsityPattern compressedSparsityPattern(1, 1);
		compressedSparsityPattern.add(0, 0);
		dealii::SparsityPattern sparsityPattern;
		sparsityPattern.copy_from(compressedSparsityPattern);
		distributed_sparse_matrix A;
		A.reinit(sparsityPattern);
		A.set(0, 0, lambda);

		DealIIWrapper<distributed_sparse_matrix, distributed_vector> tm(dt,
				name,MINRES);

		UNDISTRIBUTED_VECTOR(f,1);
		f(0) = 1;
		double t = 0;
		unsigned int n_steps = 0;
		if (7 == n){
			tmax = 0.1;
		}
		while (t < tmax) {
			if (t + dt > tmax)
				dt = tmax - t;
			t = tm.step(f, A, g, t, dt);
			if (n >= 7) { // embedded rk
				dt = tm.getDealIIEmbedded()->get_status().delta_t_guess;
			}
			++n_steps;
		}
		BOOST_CHECK_CLOSE((double ) f(0), exp(lambda * t), 1);

		cout << " " << n_steps << " steps. done." << endl;
	}
	cout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(DealIIWrapper_MultiBlock_test) {
	cout << "DealIIWrapper_MultiBlock_test..." << endl;

	for (int n = 0; n < 12; n++) {
		DealIntegratorName name = static_cast<DealIntegratorName>(n);
		cout << "  - Integrator " << n << "...";

		// solve ODE F(f) = A*f
		// with initial value f0 = [ 1 2 ]
		// and A = [ 1 -1 , 0 3 ]
		// analytic solution: f(x) = e^(A*x)
		double tmax = 1;
		double dt = 0.001;
		// build matrix
		distributed_sparse_block_matrix A;
#ifdef WITH_TRILINOS
		dealii::DynamicSparsityPattern compressedSparsityPattern(1, 1);
		compressedSparsityPattern.add(0, 0);
		A.reinit(2, 2);
		A.block(0, 0).reinit(compressedSparsityPattern);
		A.block(0, 1).reinit(A.block(0, 0));
		A.block(1, 0).reinit(A.block(0, 0));
		A.block(1, 1).reinit(A.block(0, 0));
		A.collect_sizes();
#else
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
#endif
		A.set(0, 0, 1);
		A.set(0, 1, -1);
		A.set(1, 0, 0);
		A.set(1, 1, 3);

		DealIIWrapper<distributed_sparse_block_matrix, distributed_block_vector> tm(
				dt, name,MINRES);
		// initialize block vectors
		distributed_block_vector f;
		distributed_block_vector b;
		REINIT_UNDISTRIBUTED_BLOCK_VECTOR(f, 2, 1);
		REINIT_UNDISTRIBUTED_BLOCK_VECTOR(b, 2, 1);
		f(0) = 1;
		f(1) = 2;
		b(0) = 0;
		b(1) = 0;
		double c0 = f(0);
		double c1 = f(1);
		double t = 0;

		if (7 == n){
			tmax = 0.1;
		}
		unsigned int n_steps = 0;
		while (t < tmax) {
			if (t + dt > tmax)
				dt = tmax - t;
			t = tm.step(f, A, b, t, dt);
			if (n >= 7) { // embedded rk
				dt = tm.getDealIIEmbedded()->get_status().delta_t_guess;
			}
			++n_steps;
		}

		// calculate analytic solution
		// f1(x) = c1 * exp(x) - 0.5 * c2 * exp(x) [ exp(2x) - 1]
		// f2(x) = c2 * exp(3x)
		// with c2 = f2(0) and c1 = f1(0)
		double f0 = c0 * exp(t) - 0.5 * c1 * exp(t) * (exp(2 * t) - 1.0);
		double f1 = c1 * exp(3 * t);
		// compare
		BOOST_CHECK_CLOSE((double ) f(0), f0, 1);
		BOOST_CHECK_CLOSE((double ) f(1), f1, 1);

		A.clear();

		cout << " " << n_steps << " steps. done." << endl;

	}
	cout << "done." << endl;
} /* DealIIWrapper_MultiBlock_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
