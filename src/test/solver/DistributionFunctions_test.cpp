/**
 * @file DistributionFunctions_test.cpp
 * @short
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/solver/DistributionFunctions.h"

#include "boost/test/included/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(DistributionFunctions_test)

BOOST_AUTO_TEST_CASE(DistributionFunctions_Construction_test) {
	pout << "DistributionFunctions_Construction_test..." << endl;

	/// Empty constructor
	BOOST_CHECK_NO_THROW(DistributionFunctions());

	// vectors have to be distributed, because otherwise
	// they are recognized as ghost vectors; and ghost
	// do not support writing on individual elements
	PeriodicTestDomain2D test_domain(3);
	dealii::QGaussLobatto<1> quadrature(2);
	dealii::FE_DGQArbitraryNodes<2> fe(quadrature);
	dealii::DoFHandler<2> dof_handler(*(test_domain.getMesh()));
	dof_handler.distribute_dofs(fe);

	/// Conversion-from-vector constructor
	vector<distributed_vector> f;
	for (size_t i = 0; i < 9; i++) {
		distributed_vector f_i;
		f_i.reinit((dof_handler.locally_owned_dofs()), MPI_COMM_WORLD);
		for (size_t j = 0; j < 10; j++) {
			if (f_i.in_local_range(j)) {
				f_i(j) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
						+ pow((0.5 * cos(j)), 2);
			}
		}
		f.push_back(f_i);
	}
	DistributionFunctions fD(f);
	for (size_t i = 0; i < 9; i++) {
		for (size_t j = 0; j < 10; j++) {
			if (fD.at(i).in_local_range(j)) BOOST_CHECK(f.at(i)(j) == fD.at(i)(j));
		}
	}

	// Copy constructor
	DistributionFunctions fD2(fD);
	for (size_t i = 0; i < 9; i++) {
		for (size_t j = 0; j < 10; j++) {
			if (fD.at(i).in_local_range(j)) BOOST_CHECK(fD2.at(i)(j) == fD.at(i)(j));
		}
	}

	pout << "done" << endl;
} /* DistributionFunctions_Construction_test */


/*
BOOST_AUTO_TEST_CASE(DistributionFunctions_StencilTransfer_test) {
	pout << "DistributionFunctions_StencilTransfer_test..." << endl;


	pout << "done" << endl;
} */ /* DistributionFunctions_StencilTransfer_test */

BOOST_AUTO_TEST_SUITE_END()

