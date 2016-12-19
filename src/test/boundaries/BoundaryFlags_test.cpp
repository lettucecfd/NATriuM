/**
 * @file BoundaryFlags_test.cpp
 * @short Unit tests for all functions ins BoundaryTools.h
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/boundaries/BoundaryFlags.h"

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(BoundaryFlags_test)

BOOST_AUTO_TEST_CASE(BoundaryFlags_Operators_test) {

	pout << "BoundaryFlags_Operators_test..." << endl;

	BoundaryFlags flags = only_distributions;
	BOOST_CHECK_EQUAL(flags, only_distributions);

	flags |= boundary_rho;
	BOOST_CHECK_EQUAL(flags, boundary_rho);

	flags |= boundary_drho_dt;
	BOOST_CHECK(boundary_rho & flags);
	BOOST_CHECK(boundary_drho_dt & flags);
	BOOST_CHECK(not (boundary_u & flags));

	pout << "done." << endl;

} /* BoundaryFlags_Operators_test */

BOOST_AUTO_TEST_CASE(BoundaryFlags_PrescribedQuantities_test) {

	pout << "BoundaryFlags_PrescribedQuantities_test..." << endl;

	dealii::Tensor<1,2> u;
	dealii::Point<2> p;
	PrescribedQuantities<2> pu(u);
	BOOST_CHECK_EQUAL(pu.getPrescribedValues(), boundary_u);
	BOOST_CHECK(not pu.getPressure());
	BOOST_CHECK_EQUAL( pu.getVelocity()->value(p, 1), 0.0);

	PrescribedQuantities<2> pp (1.0);
	BOOST_CHECK(pp.getPrescribedValues() == boundary_p);
	BOOST_CHECK_EQUAL(pp.getPressure()->value(p), 1.0);
	BOOST_CHECK(not pp.getVelocity());

	pout << "done." << endl;

} /* BoundaryFlags_PrescribedQuantities_test */

BOOST_AUTO_TEST_SUITE_END()
