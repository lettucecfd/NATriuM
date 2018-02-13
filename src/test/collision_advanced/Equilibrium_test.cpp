/*
 * Equilibrium_test.cpp
 *
 *  Created on: 20.01.2017
 *      Author: natrium
 */
#include "natrium/stencils/D2Q9.h"
#include "natrium/collision/BGKStandard.h"

#include <math.h>
#include <array>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/collision_advanced/Equilibria.h"
#include "natrium/collision_advanced/AuxiliaryCollisionFunctions.h"
#include "natrium/solver/SolverConfiguration.h"

#include "deal.II/fe/fe_dgq.h"

#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/benchmarks/PeriodicTestDomain3D.h"
#include "../problemdescription/TaylorGreenTest2D.h"

using namespace natrium;


BOOST_AUTO_TEST_SUITE(Equilibrium_test_suite)

BOOST_AUTO_TEST_CASE(Equilibrium_test) {



	// initialize distributions with arbitrary components
//vector<distributed_vector> f;
double rho = 1.0;
		//rho.compress(dealii::VectorOperation::add);
	vector<distributed_vector> u;
/*	for (size_t i = 0; i < 9; i++) {
		distributed_vector f_i;
				f_i(0) = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
						+ pow((0.5 * cos(0.5)), 2);

		//f_i.compress(dealii::VectorOperation::add);
		f.push_back(f_i);
	}*/

	//DistributionFunctions f_new(f);

	double velocities[2]={0.1,0.2};

	BGKEquilibrium<2,9> eq;
	std::array<double,9> feq;
	D2Q9 d2q9(1.0);
	double cs2=1./3.;
	double scaling = 1.0;
	double dt = 0.1;

	double viscosity = 1.0;


	SolverConfiguration cfg;
	TaylorGreenTest2D tgv(0.1,1);
	GeneralCollisionData<2,9> prams(cfg, tgv, scaling, viscosity, d2q9,
			cs2 , dt);

	prams.velocity[0]=velocities[0];
	prams.velocity[1]=velocities[1];
	prams.density = rho;


	eq.calc(feq,prams);

	/*for (int i =0;i<9;i++)
	{
		cout << feq[i] << endl;
	}*/

	//cout << "Dichte v:" << prams.density << endl;
	//cout << "Geschwindigkeit vorher:" << prams.velocity[0] << endl;
	//cout << "Geschwindigkeit vorher:" << prams.velocity[1] << endl;



	BOOST_CHECK(feq[1]-feq[3]+feq[5]-feq[6]-feq[7]+feq[8]-prams.velocity[0]<10e-6);
	BOOST_CHECK(feq[1]-feq[3]+feq[5]-feq[6]-feq[7]+feq[8]-prams.velocity[0]<10e-6);





} // test_case

BOOST_AUTO_TEST_SUITE_END()




