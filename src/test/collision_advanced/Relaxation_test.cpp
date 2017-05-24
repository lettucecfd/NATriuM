


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

BOOST_AUTO_TEST_CASE(Relaxation_test) {



	// initialize distributions with arbitrary components
	double test_f[9];

double rho = 1.0;
		//rho.compress(dealii::VectorOperation::add);
	vector<distributed_vector> u;
	for (size_t i = 0; i < 9; i++) {
				test_f[i] = 1.5 + sin(1.5 * i) + 0.001 + i / (i + 1)
						+ pow((0.5 * cos(0.5)), 2);
	}

	double cs2=1./3.;
	double scaling = 1.0;
	double dt = 0.1;

	double viscosity = 1.0;
	
	D2Q9 d2q9(1.0);
	GeneralCollisionData<2,9> prams(scaling, viscosity, d2q9,
			cs2 , dt);


	cout << "Tau: " << calculateTauFromNu(viscosity,cs2,dt);


	Regularized<2, 9, BGKEquilibrium>::SpecificCollisionData data(prams);
	prams.density = calculateDensity<9>(test_f);
	calculateVelocity<2,9>(test_f,prams.velocity,1.0,prams.density, prams);



	Regularized<2,9,BGKEquilibrium> test;



	test.relax(test_f,prams,data);
	cout << "Dichte v:" << prams.density << endl;
	cout << "Geschwindigkeit vorher:" << prams.velocity[0] << endl;
	cout << "Geschwindigkeit vorher:" << prams.velocity[1] << endl;

	calculateVelocity<2,9>(test_f,prams.velocity,1.0,prams.density, prams);
	cout << "Dichte n:" << prams.density << endl;
	cout << "Geschwindigkeit nachher:" << prams.velocity[0] << endl;
	cout << "Geschwindigkeit nachher:" << prams.velocity[1] << endl;

	for (int i=0;i<9;i++)
	{
		cout << prams.weight[i] << endl;
	}

	cout << (test_f[1]-test_f[3]+test_f[5]-test_f[6]-test_f[7]+test_f[8]) / prams.density << endl << prams.density << endl;

	BOOST_CHECK((test_f[1]-test_f[3]+test_f[5]-test_f[6]-test_f[7]+test_f[8])/prams.density-prams.velocity[0]<10e-6);








	} // test_case

BOOST_AUTO_TEST_SUITE_END()


