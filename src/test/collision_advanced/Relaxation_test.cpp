


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

#include "deal.II/fe/fe_dgq.h"

#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/benchmarks/PeriodicTestDomain3D.h"

namespace natrium{


//BOOST_AUTO_TEST_SUITE(Relaxation_test_suite)

BOOST_AUTO_TEST_CASE(Relaxation_test_suite) {



	// initialize distributions with arbitrary components
	double test_f[9];

double rho = 1.0;
		//rho.compress(dealii::VectorOperation::add);
	vector<distributed_vector> u;
	for (size_t i = 0; i < 9; i++) {
				test_f[i] = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
						+ pow((0.5 * cos(0.5)), 2);
	}

	double velocities[2];
	double density;
	density = calculateDensity<9>(test_f);
	calculateVelocity<2,9>(test_f,velocities,1.0,density);

	CollisionParameters<2,9> prams;
	prams.cs2=1./3.;
	prams.scaling = 1.0;
	prams.velocity[0]=velocities[0];
	prams.velocity[1]=velocities[1];
	prams.density = density;
	prams.tau = 1.0;

	BGKCollision<2,9,BGKEquilibrium> test;



	test.relax(test_f,prams);

	cout << (test_f[1]-test_f[3]+test_f[5]-test_f[6]-test_f[7]+test_f[8]) / density << endl << density << endl;

	BOOST_CHECK((test_f[1]-test_f[3]+test_f[5]-test_f[6]-test_f[7]+test_f[8])/density-velocities[0]<10e-6);








	} // test_case

BOOST_AUTO_TEST_SUITE_END()


