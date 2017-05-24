#include "natrium/collision_advanced/AuxiliaryMRTFunctions.h"
#include "boost/test/unit_test.hpp"
#include "natrium/utilities/BasicNames.h"
#include "natrium/solver/SolverConfiguration.h"

#include <array>

using std::array;
using namespace natrium;
using namespace natrium::AuxiliaryMRTFunctions;

BOOST_AUTO_TEST_SUITE(AuxiliaryMRTFunctions_test)

BOOST_AUTO_TEST_CASE(AuxiliaryMRTFunctions_Invert_test) {

	pout << "AuxiliaryMRTFunctions_Invert_test..." << endl;
	array<array<double, 9>, 9> result9 = { };
	array<array<double, 19>, 19> result19 = { };

	// Dellar D2Q9
	matrix_matrix_product(MRTDellarD2Q9::moment_trafo,
			MRTDellarD2Q9::inverse_trafo, result9);
	for (size_t i = 0; i < 9; i++) {
		for (size_t j = 0; j < 9; j++) {
			BOOST_CHECK_SMALL(result9[i][j] - double(i == j), 1e-12);
		}
	}

	// Lallemand D2Q9
	matrix_matrix_product(MRTLallemandD2Q9::moment_trafo,
			MRTLallemandD2Q9::inverse_trafo, result9);
	for (size_t i = 0; i < 9; i++) {
		for (size_t j = 0; j < 9; j++) {
			BOOST_CHECK_SMALL(result9[i][j] - double(i == j), 1e-12);
		}
	}

	// DHumieres D3Q19
	matrix_matrix_product(MRTDHumieresD3Q19::moment_trafo,
			MRTDHumieresD3Q19::inverse_trafo, result19);
	for (size_t i = 0; i < 9; i++) {
		for (size_t j = 0; j < 9; j++) {
			BOOST_CHECK_SMALL(result19[i][j] - double(i == j), 1e-12);
		}
	}

	pout << "done." << endl;
} /* AuxiliaryMRTFunctions_Invert_test */



BOOST_AUTO_TEST_SUITE_END()

