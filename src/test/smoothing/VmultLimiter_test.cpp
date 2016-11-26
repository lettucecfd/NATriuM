/*
 * VmultLimiter_test.cpp
 *
 *  Created on: 29.08.2016
 *      Author: akraem3m
 */

#include "natrium/smoothing/VmultLimiter.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/base/index_set.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/advection/SemiLagrangian.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(VmultLimiter_test)

BOOST_AUTO_TEST_CASE(VmultLimiter_apply_test) {

	pout << "VmultLimiter_apply_test..." << endl;


	// make a matrix
	PeriodicTestDomain2D domain(3);
	domain.refineAndTransform();
	double dt = 0.03/(1.0*0.125); // so that dt xi dx = 0.03; (should be unstable)

	SemiLagrangian<2> semi(domain, 3, boost::make_shared<D2Q9>(), dt);
	semi.setupDoFs();
	//semi.setDeltaT(dt);
	semi.reassemble();
	const dealii::TrilinosWrappers::SparseMatrix& M = semi.getSystemMatrix().block(0,0);
	dealii::TrilinosWrappers::MPI::Vector s(M.locally_owned_domain_indices());
	dealii::TrilinosWrappers::MPI::Vector t(M.locally_owned_range_indices());


	// fill source with zeros, where Mij < 0 and ones everywhere else
	int i = *(M.locally_owned_range_indices().begin());
	s = 1;
	for (size_t j = 0; j < M.n(); j++){
		if (M.el(i,j) < 0)
			s(j) = 0;
	}
	//s.compress();

	// make sure target(index_0) > 1
	M.vmult(t, s);
	double ti = t(i);
	BOOST_CHECK_GT(ti, 1);
	// apply limiter
	VmultLimiter::apply(M,t,s);

	// make sure target(index_0) == 1
	ti = t(i);
	BOOST_CHECK_LE(ti, 1);


	pout << "done" << endl;
} /*VmultLimiter_apply_test*/

BOOST_AUTO_TEST_SUITE_END()

