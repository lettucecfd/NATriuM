/**
 * @file DistributionFunctions_test.cpp
 * @short
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/solver/DistributionFunctions.h"

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(DistributionFunctions_test)

BOOST_AUTO_TEST_CASE(DistributionFunctions_Construction_test) {
	pout << "DistributionFunctions_Construction_test..." << endl;

	/// Empty constructor
	BOOST_CHECK_NO_THROW(DistributionFunctions());

	/// Conversion-from-vector constructor
	vector<distributed_vector> f;
	for (size_t i = 0; i < 9; i++){
		UNDISTRIBUTED_VECTOR(f_i,10);
		for (size_t j = 0; j < 10; j++){
			f_i(j) = 1.5 + sin(1.5*i)+0.001+i/(i+1) + pow((0.5*cos(j)),2);
		}
		f.push_back(f_i);
	}
	DistributionFunctions fD(f);
	for (size_t i = 0; i < 9; i++){
		for (size_t j = 0; j < 10; j++){
			BOOST_CHECK(f.at(i)(j) == fD.at(i)(j));
		}
	}

	// Copy constructor
	DistributionFunctions fD2(fD);
	for (size_t i = 0; i < 9; i++){
		for (size_t j = 0; j < 10; j++){
			BOOST_CHECK(fD2.at(i)(j) == fD.at(i)(j));
		}
	}


	pout << "done" << endl;
} /* DistributionFunctions_Construction_test */


BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
