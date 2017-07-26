/*
 * EntropicStabilized_test.cpp
 *
 *  Created on: 11.05.2017
 *      Author: akraem3m
 */

#include "natrium/collision/EntropicStabilized.h"

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include <array>

using std::array;

namespace natrium {

using namespace natrium;

BOOST_AUTO_TEST_SUITE(EntropicStabilized_test)

BOOST_AUTO_TEST_CASE(EntropicStabilizedDensityMomentum_test) {
	pout << "EntropicStabilizedDensityMomentum_test..." << endl;

	D2Q9 st;
	array<double, 9> f = {{ 1.0, 0.5, 0.1, 0.2, 0.5, 0.2, 0.1, 1.0, 0.5 }};
	array<double, 2> j = { };
	double rho = density<9>(f);
	momentum<2, 9>(f, j, st.getDirections());

	BOOST_CHECK_CLOSE(rho, 4.1, 1e-10);
	BOOST_CHECK_CLOSE(j[0], -0.1, 1e-10);
	BOOST_CHECK_CLOSE(j[1], -1.6, 1e-10);

	// with scaled stencil
	D2Q9 st2(2);
	array<double, 2> j2 = { };
	double rho2 = density<9>(f);
	momentum<2, 9>(f, j2, st2.getDirections());

	BOOST_CHECK_CLOSE(rho2, 4.1, 1e-10);
	BOOST_CHECK_CLOSE(j2[0], 2*-0.1, 1e-10);
	BOOST_CHECK_CLOSE(j2[1], 2*-1.6, 1e-10);
	pout << "done" << endl;
} /* EntropicStabilizedDensityMomentum_test */

BOOST_AUTO_TEST_CASE(EntropicStabilized_collideBGK_test) {
	pout << "EntropicStabilized_collideBGK_test..." << endl;

	D2Q9 st;
	array<double, 9> f = {{ 1.0, 0.5, 0.1, 0.2, 0.5, 0.2, 0.1, 1.0, 0.5 }};
	array<double, 9> feq = { };
	array<double, 2> j = { };
	array<double, 2> u = { };
	double rho = density<9>(f);
	double tau = 1.0;
	momentum<2, 9>(f, j, st.getDirections());
	u[0] = j[0] / rho;
	u[1] = j[1] / rho;

	// collide
	collideBGK<2, 9>(f, rho, u, tau, st, feq);

	// check conservation
	array<double, 2> j_out = { };
	double rho_out = density<9>(f);
	momentum<2, 9>(f, j_out, st.getDirections());
	BOOST_CHECK_CLOSE(rho_out, rho, 1e-10);
	BOOST_CHECK_CLOSE(j_out[0], j[0], 1e-10);
	BOOST_CHECK_CLOSE(j_out[1], j[1], 1e-10);

	// check full relaxation
	for (size_t i = 0; i < 9; i++) {
		BOOST_CHECK_CLOSE(f[i], feq[i], 1e-10);
	}

	// check overrelaxation
	double f0 = f[0];
	collideBGK<2, 9>(f, rho, u, tau, st, feq);
	BOOST_CHECK_CLOSE(f[0], f0 + 1. / tau * (f0 - feq[0]), 1e-10);

	pout << "done" << endl;
} /* EntropicStabilized_collideBGK_test */

BOOST_AUTO_TEST_CASE(EntropicStabilizedCollideOneD2Q9_test) {
	pout << "EntropicStabilizedCollideOneD2Q9_test..." << endl;

	boost::shared_ptr<Stencil> st = boost::make_shared<D2Q9>();
	double tau = 1.0;
	double dt = 1.0;
	EStCollisionData<2, 9> _(*st, tau);
	_.f_i = {{1.0, 0.5, 0.1, 0.2, 0.5, 0.2, 0.1, 1.0, 0.5}};

	// calc density and momentum (as input to the collision function)
	_.rho_i = density<9>(_.f_i);
	momentum<2, 9>(_.f_i, _.j_i, st->getDirections());
	_.u_i[0] = _.j_i[0] / _.rho_i;
	_.u_i[1] = _.j_i[1] / _.rho_i;

	// collide
	EntropicStabilized est(tau - 0.5, dt, st);
	est.collideOne<2, 9>(_);

	// calculate density and momentum
	array<double, 2> j = { };
	double rho = density<9>(_.f_i);
	momentum<2, 9>(_.f_i, j, st->getDirections());

	BOOST_CHECK_CLOSE(rho, _.rho_i, 1e-10);
	BOOST_CHECK_CLOSE(j[0], _.j_i[0], 1e-10);
	BOOST_CHECK_CLOSE(j[1], _.j_i[1], 1e-10);

	pout << "done" << endl;
} /* EntropicStabilizedCollideOneD2Q9_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
