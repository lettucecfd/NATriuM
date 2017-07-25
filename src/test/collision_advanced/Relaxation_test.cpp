#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/BGKStandard.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/collision_advanced/Equilibria.h"
#include "natrium/collision_advanced/CollisionSchemes.h"
#include "natrium/collision_advanced/AuxiliaryCollisionFunctions.h"
#include "natrium/solver/SolverConfiguration.h"

#include "deal.II/fe/fe_dgq.h"

#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/benchmarks/PeriodicTestDomain3D.h"
#include "../problemdescription/TaylorGreenTest2D.h"
#include "natrium/benchmarks/PoiseuilleFlow3D.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(Relaxation_test_suite)

BOOST_AUTO_TEST_CASE(Relaxation_Regularized_test) {

	pout << "Relaxation_Regularized_test..." << endl;

	// initialize distributions with arbitrary components
	std::array<double, 9> test_f;

	double rho = 1.0;
	for (size_t i = 0; i < 9; i++) {
		test_f[i] = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
				+ pow((0.5 * cos(0.5)), 2);
	}

	double cs2 = 1. / 3.;
	double scaling = 1.0;
	double dt = 0.1;

	double viscosity = 1.0;

	D2Q9 d2q9(1.0);
	SolverConfiguration cfg;
	TaylorGreenTest2D tgv(0.1,1);
	GeneralCollisionData<2, 9> prams(cfg, tgv, scaling, viscosity, d2q9, cs2, dt);

	Regularized<2, 9, BGKEquilibrium>::SpecificCollisionData data(prams);

	prams.density = calculateDensity<9>(test_f);
	calculateVelocity<2, 9>(test_f, prams.velocity, prams.density, prams);

	Regularized<2, 9, BGKEquilibrium> test;

	test.relax(test_f, prams, data);

	calculateVelocity<2, 9>(test_f, prams.velocity,  prams.density, prams);

	BOOST_CHECK(
			(test_f[1] - test_f[3] + test_f[5] - test_f[6] - test_f[7]
					+ test_f[8]) / prams.density - prams.velocity[0] < 10e-6);

	pout << "done." << endl;

} /* Relaxation_Regularized_test */

// ==============================================================
// ======================== MRT =================================
// ==============================================================

BOOST_AUTO_TEST_CASE(Relaxation_MRT_D2Q9_test) {
	// test conservation of mass and momentum
	pout << "Relaxation_MRT_D2Q9_test..." << endl;

	// initialize distributions with arbitrary components
	std::array<double, 9> test_f;
	double rho = 1.0;
	std::array<double, 2> u = { };
	for (size_t i = 0; i < 9; i++) {
		test_f[i] = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
				+ pow((0.5 * cos(0.5)), 2);
	}
	double scaling = 4.0;
	double cs2 = scaling * scaling / 3.;
	double dt = 0.1;
	double viscosity = 1.0;

	D2Q9 d2q9(scaling);
	SolverConfiguration cfg;

	// ======================== Dellar D2Q9 =====================
	cfg.setMRTBasis(DELLAR_D2Q9);
	TaylorGreenTest2D tgv(0.1,1);
	GeneralCollisionData<2, 9> prams(cfg, tgv, scaling, viscosity, d2q9, cs2, dt);

	MultipleRelaxationTime<2, 9, BGKEquilibrium>::SpecificCollisionData data(
			prams);
	prams.density = calculateDensity<9>(test_f);
	calculateVelocity<2, 9>(test_f, prams.velocity,  prams.density,
			prams);
	rho = prams.density;
	u[0] = prams.velocity[0];
	u[1] = prams.velocity[1];

	MultipleRelaxationTime<2, 9, BGKEquilibrium> test;
	test.relax(test_f, prams, data);
	prams.density = calculateDensity<9>(test_f);
	calculateVelocity<2, 9>(test_f, prams.velocity,  prams.density,
			prams);

	BOOST_CHECK_CLOSE(rho, prams.density, 1e-10);
	BOOST_CHECK_CLOSE(u[0], prams.velocity[0], 1e-10);
	BOOST_CHECK_CLOSE(u[1], prams.velocity[1], 1e-10);

	// ======================== Lallemand D2Q9 =====================
	cfg.setMRTBasis(DELLAR_D2Q9);
	GeneralCollisionData<2, 9> gendata(cfg, tgv, scaling, viscosity, d2q9, cs2, dt);

	MultipleRelaxationTime<2, 9, BGKEquilibrium>::SpecificCollisionData data2(
			gendata);
	gendata.density = calculateDensity<9>(test_f);
	calculateVelocity<2, 9>(test_f, gendata.velocity,  gendata.density,
			gendata);
	rho = gendata.density;
	u[0] = gendata.velocity[0];
	u[1] = gendata.velocity[1];

	MultipleRelaxationTime<2, 9, BGKEquilibrium> test2;
	test2.relax(test_f, gendata, data2);
	gendata.density = calculateDensity<9>(test_f);
	calculateVelocity<2, 9>(test_f, gendata.velocity,  gendata.density,
			gendata);

	BOOST_CHECK_CLOSE(rho, gendata.density, 1e-10);
	BOOST_CHECK_CLOSE(u[0], gendata.velocity[0], 1e-10);
	BOOST_CHECK_CLOSE(u[1], gendata.velocity[1], 1e-10);

	pout << "done." << endl;

} /* Relaxation_MRT_D2Q9_test */

BOOST_AUTO_TEST_CASE(Relaxation_MRT_D3Q19_test) {
	// test conservation of mass and momentum
	pout << "Relaxation_MRT_D3Q19_test..." << endl;

	// initialize distributions with arbitrary components
	std::array<double, 19> test_f;
	double rho = 1.0;
	std::array<double, 3> u = { };
	for (size_t i = 0; i < 19; i++) {
		test_f[i] = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
				+ pow((0.5 * cos(0.5)), 2);
	}
	double scaling = 4.0;
	double cs2 = scaling * scaling / 3.;
	double dt = 0.1;
	double viscosity = 1.0;

	D3Q19 d3q19(scaling);
	SolverConfiguration cfg;

	// ======================== DHumieres D3Q19 =====================
	cfg.setMRTBasis(DHUMIERES_D3Q19);
	PoiseuilleFlow3D pois(0.1,1); // just a dummy
	GeneralCollisionData<3, 19> prams(cfg, pois, scaling, viscosity, d3q19, cs2, dt);

	MultipleRelaxationTime<3, 19, BGKEquilibrium>::SpecificCollisionData data(
			prams);
	prams.density = calculateDensity<19>(test_f);
	calculateVelocity<3, 19>(test_f, prams.velocity,  prams.density,
			prams);
	rho = prams.density;
	u[0] = prams.velocity[0];
	u[1] = prams.velocity[1];
	u[2] = prams.velocity[2];

	MultipleRelaxationTime<3, 19, BGKEquilibrium> test;
	test.relax(test_f, prams, data);
	prams.density = calculateDensity<19>(test_f);
	calculateVelocity<3, 19>(test_f, prams.velocity,  prams.density,
			prams);

	BOOST_CHECK_CLOSE(rho, prams.density, 1e-10);
	BOOST_CHECK_CLOSE(u[0], prams.velocity[0], 1e-10);
	BOOST_CHECK_CLOSE(u[1], prams.velocity[1], 1e-10);
	BOOST_CHECK_CLOSE(u[2], prams.velocity[2], 1e-10);

	pout << "done." << endl;

} /* Relaxation_MRT_D3Q19_test */



BOOST_AUTO_TEST_CASE(Relaxation_Equiv_MRT_Reg_test) {
	// test conservation of mass and momentum
	// Dellar D2Q9 with full relaxation should be equivalent to Regularized LB

	pout << "Relaxation_Equiv_MRT_Reg_test..." << endl;

	// initialize distributions with arbitrary components
	std::array<double, 9> test_f;
	std::array<double, 9> test_f2;
	double rho = 1.0;
	std::array<double, 2> u = { };
	for (size_t i = 0; i < 9; i++) {
		test_f[i] = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
				+ pow((0.5 * cos(0.5)), 2);
		test_f2[i] = test_f[i];
	}
	double scaling = 4.0;
	double cs2 = scaling * scaling / 3.;
	double dt = 0.1;
	double viscosity = 1.0;

	D2Q9 d2q9(scaling);
	SolverConfiguration cfg;
	TaylorGreenTest2D tgv(0.1,1); // just a dummy

	// ======================== Dellar D2Q9 =====================
	cfg.setMRTBasis(DELLAR_D2Q9);
	GeneralCollisionData<2, 9> prams(cfg, tgv, scaling, viscosity, d2q9, cs2, dt);

	MultipleRelaxationTime<2, 9, BGKEquilibrium>::SpecificCollisionData data(
			prams);
	prams.density = calculateDensity<9>(test_f);
	calculateVelocity<2, 9>(test_f, prams.velocity,  prams.density,
			prams);

	MultipleRelaxationTime<2, 9, BGKEquilibrium> test;
	test.relax(test_f, prams, data);

	// ======================== Regularized =====================

	GeneralCollisionData<2, 9> prams2(cfg, tgv, scaling, viscosity, d2q9, cs2, dt);
	Regularized<2, 9, BGKEquilibrium>::SpecificCollisionData data2(
			prams2);
	prams2.density = calculateDensity<9>(test_f2);
	calculateVelocity<2, 9>(test_f2, prams2.velocity,  prams2.density,
			prams2);

	Regularized<2, 9, BGKEquilibrium> test2;
	test2.relax(test_f2, prams2, data2);

	for (size_t i = 0; i < 9; i++){
		BOOST_CHECK_CLOSE(test_f2[i], test_f[i], 1e-10);
	}

	pout << "done." << endl;
} /* Relaxation_MRT_D3Q19_test */

BOOST_AUTO_TEST_SUITE_END()

