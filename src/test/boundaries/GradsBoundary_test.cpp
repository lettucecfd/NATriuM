#include "natrium/boundaries/GradsBoundary.h"

#include <string>

#include "boost/test/unit_test.hpp"

#include "deal.II/base/tensor.h"

#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(GradsBoundary_test)

BOOST_AUTO_TEST_CASE(GradsBoundary_Constructor_test) {

	pout << "GradsFunction_Call_test..." << endl;

	dealii::Tensor<2, 2> P2;
	dealii::Tensor<2, 3> P3;
	dealii::Tensor<1, 2> j2;
	dealii::Tensor<1, 3> j3;
	double rho = 1;
	vector<double> f;

	f.resize(9);
	GradsFunction<2>(f, D2Q9(), rho, j2, P2);
	f.resize(15);
	GradsFunction<3>(f, D3Q15(), rho, j3, P3);
	f.resize(19);
	GradsFunction<3>(f, D3Q19(), rho, j3, P3);
	f.resize(27);
	GradsFunction<3>(f, D3Q27(), rho, j3, P3);

	pout << "done." << endl;

} /* GradsFunction_Construction_test */

BOOST_AUTO_TEST_CASE(GradsBoundary_Velocity_test) {

} /* GradsFunction_Moments_test */

BOOST_AUTO_TEST_SUITE_END()
} /* namespace natrium */

