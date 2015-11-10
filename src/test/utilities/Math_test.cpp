/**
 * @file Math_test.cpp
 * @short
 * @date 23.05.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/Math.h"

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"

namespace natrium{

BOOST_AUTO_TEST_SUITE(Math_test)

BOOST_AUTO_TEST_CASE(Math_VelocityNorm_test){
	pout << "Math_VelocityNorm_test..." << endl;

	// 1d case -> failure
	vector<distributed_vector> v;
	UNDISTRIBUTED_VECTOR(vx,2);
	vx(0) = 0.3;
	vx(1) = 0.1;
	v.push_back(vx);
	//BOOST_CHECK_THROW(Math::maxVelocityNorm(v), std::exception);
	//BOOST_CHECK_THROW(Math::velocity2Norm(v), std::exception);

	// 2d case
	UNDISTRIBUTED_VECTOR(vy,2);
	vy(0) = 0.4;
	vy(1) = 0.2;
	v.push_back(vy);
	BOOST_CHECK_CLOSE(Math::maxVelocityNorm(v), 0.5, 1e-15);
	BOOST_CHECK_CLOSE(Math::velocity2Norm(v), sqrt(30)/10., 1e-15);

	// 3d case
	UNDISTRIBUTED_VECTOR(vz,2);
	vz(0) = 0.2;
	vz(1) = - sqrt(76)/10;
	v.push_back(vz);
	BOOST_CHECK_CLOSE(Math::maxVelocityNorm(v), 0.9, 1e-13);
	BOOST_CHECK_CLOSE(Math::velocity2Norm(v), sqrt(110)/10., 1e-13);

	// 4d case -> failure
	UNDISTRIBUTED_VECTOR(va,2);
	v.push_back(va);
	//BOOST_CHECK_THROW(Math::maxVelocityNorm(v), std::exception);
	//BOOST_CHECK_THROW(Math::velocity2Norm(v), std::exception);

	pout << "done" <<endl;
}

BOOST_AUTO_TEST_SUITE_END()

}/*namespace natrium*/
