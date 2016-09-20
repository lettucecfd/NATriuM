#include "natrium/boundaries/GradsFunction.h"

#include <string>

#include "boost/test/unit_test.hpp"

#include "deal.II/base/tensor.h"

#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {
namespace BoundaryTools {

BOOST_AUTO_TEST_SUITE(GradsFunction_test)

BOOST_AUTO_TEST_CASE(GradsFunction_Call_test) {

	pout << "GradsFunction_Call_test..." << endl;

	dealii::Tensor<2, 2> P2;
	dealii::Tensor<2, 3> P3;
	dealii::Tensor<1, 2> j2;
	dealii::Tensor<1, 3> j3;
	double rho = 1;
	vector<double> f;

	f.resize(9);
	GradsFunction<D2Q9, 2>(f, D2Q9(), rho, j2, P2);
	f.resize(15);
	GradsFunction<D3Q15, 3>(f, D3Q15(), rho, j3, P3);
	f.resize(19);
	GradsFunction<D3Q19, 3>(f, D3Q19(), rho, j3, P3);
	f.resize(27);
	GradsFunction<D3Q27, 3>(f, D3Q27(), rho, j3, P3);

	pout << "done." << endl;

} /* GradsFunction_Construction_test */

BOOST_AUTO_TEST_CASE(GradsFunction_Moments_test) {

	pout << "GradsFunction_Moments_test..." << endl;

	dealii::Tensor<2, 2> P2;
	P2[0][0] = 1.5;
	P2[0][1] = 1.1;
	P2[1][0] = 1.1;
	P2[1][1] = 0.3;
	dealii::Tensor<2, 3> P3;
	P3[0][0] = 1.5;
	P3[0][1] = 1.1;
	P3[0][2] = 0.9;
	P3[1][0] = 1.1;
	P3[1][1] = 1.4;
	P3[1][2] = 1.3;
	P3[2][0] = 0.9;
	P3[2][1] = 1.3;
	P3[2][2] = 1.2;
	dealii::Tensor<1, 2> j2;
	j2[0] = 1;
	j2[1] = 0.31;
	dealii::Tensor<1, 3> j3;
	j3[0] = 1;
	j3[1] = 0.31;
	j3[2] = 0.41;
	double rho = 0.9;
	vector<double> f;

	double rho_out;
	double j_out;
	double P_out;

	// D2Q9
	f.resize(9);
	D2Q9 dq(12);
	GradsFunction<D2Q9, 2>(f, dq, rho, j2, P2);
	rho_out = 0;
	for (size_t alpha = 0; alpha < f.size(); alpha++) {
		rho_out += f.at(alpha);
	}
	BOOST_CHECK_CLOSE(rho_out, rho, 1e-10);
	for (size_t i = 0; i < 2; i++) {
		j_out = 0;
		for (size_t alpha = 0; alpha < f.size(); alpha++) {
			j_out += f.at(alpha) * dq.getDirection(alpha)[i];
		}
		BOOST_CHECK_CLOSE(j_out, j2[i], 1e-10);
	}
	for (size_t i = 0; i < 2; i++) {
		for (size_t j = 0; j < 2; j++) {
			P_out = 0;
			for (size_t alpha = 0; alpha < f.size(); alpha++) {
				P_out += f.at(alpha) * dq.getDirection(alpha)[i] * dq.getDirection(alpha)[j];
			}
			BOOST_CHECK_CLOSE(P_out, P2[i][j], 1e-10);
		}
	}

	// D3Q15
	f.resize(15);
	D3Q15 d3q15(1.9);
	GradsFunction<D3Q15, 3>(f, d3q15, rho, j3, P3);
	rho_out = 0;
	for (size_t alpha = 0; alpha < f.size(); alpha++) {
		rho_out += f.at(alpha);
	}
	BOOST_CHECK_CLOSE(rho_out, rho, 1e-10);
	for (size_t i = 0; i < 3; i++) {
		j_out = 0;
		for (size_t alpha = 0; alpha < f.size(); alpha++) {
			j_out += f.at(alpha) * d3q15.getDirection(alpha)[i];
		}
		BOOST_CHECK_CLOSE(j_out, j3[i], 1e-10);
	}
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			P_out = 0;
			for (size_t alpha = 0; alpha < f.size(); alpha++) {
				P_out += f.at(alpha) * d3q15.getDirection(alpha)[i] * d3q15.getDirection(alpha)[j];
			}
			BOOST_CHECK_CLOSE(P_out, P3[i][j], 1e-10);
		}
	}

	// D3Q19
	f.resize(19);
	D3Q19 d3q19(0.012);
	GradsFunction<D3Q19, 3>(f, d3q19, rho, j3, P3);
	rho_out = 0;
	for (size_t alpha = 0; alpha < f.size(); alpha++) {
		rho_out += f.at(alpha);
	}
	BOOST_CHECK_CLOSE(rho_out, rho, 1e-9);
	for (size_t i = 0; i < 3; i++) {
		j_out = 0;
		for (size_t alpha = 0; alpha < f.size(); alpha++) {
			j_out += f.at(alpha) * d3q19.getDirection(alpha)[i];
		}
		BOOST_CHECK_CLOSE(j_out, j3[i], 1e-10);
	}
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			P_out = 0;
			for (size_t alpha = 0; alpha < f.size(); alpha++) {
				P_out += f.at(alpha) * d3q19.getDirection(alpha)[i] * d3q19.getDirection(alpha)[j];
			}
			BOOST_CHECK_CLOSE(P_out, P3[i][j], 1e-10);
		}
	}

	// D3Q27
	f.resize(27);
	D3Q27 d3q27(5.2);
	GradsFunction<D3Q27, 3>(f, d3q27, rho, j3, P3);
	rho_out = 0;
	for (size_t alpha = 0; alpha < f.size(); alpha++) {
		rho_out += f.at(alpha);
	}
	BOOST_CHECK_CLOSE(rho_out, rho, 1e-10);
	for (size_t i = 0; i < 3; i++) {
		j_out = 0;
		for (size_t alpha = 0; alpha < f.size(); alpha++) {
			j_out += f.at(alpha) * d3q27.getDirection(alpha)[i];
		}
		BOOST_CHECK_CLOSE(j_out, j3[i], 1e-10);
	}
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			P_out = 0;
			for (size_t alpha = 0; alpha < f.size(); alpha++) {
				P_out += f.at(alpha) * d3q27.getDirection(alpha)[i] * d3q27.getDirection(alpha)[j];
			}
			BOOST_CHECK_CLOSE(P_out, P3[i][j], 1e-10);
		}
	}

	pout << "done." << endl;

} /* GradsFunction_Moments_test */

BOOST_AUTO_TEST_SUITE_END()} /* namepace BoundaryTools */
} /* namespace natrium */

