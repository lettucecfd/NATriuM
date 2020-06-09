/*
 * Equilibrium_test.cpp
 *
 *  Created on: 20.01.2017
 *      Author: natrium
 */
#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D2Q25H.h"
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
#include "natrium/stencils/D3Q77.h"

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



        std::array<double,2> v_post;
        calculateVelocity<2,9>(feq,v_post,calculateDensity<9>(feq),prams);

        BOOST_CHECK_CLOSE(v_post[0],prams.velocity[0],10e-6);
        BOOST_CHECK_CLOSE(v_post[1],prams.velocity[1],10e-6);






} // test_case

    BOOST_AUTO_TEST_CASE(Equilibrium_D2Q25H_test) {



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

        QuarticEquilibrium<2,25> eq;
        std::array<double,25> feq;
        D2Q25H d2q25h(1.0);
        double cs2=1./3.;
        double scaling = 1.0;
        double dt = 0.1;

        double viscosity = 1.0;


        SolverConfiguration cfg;
        TaylorGreenTest2D tgv(0.1,1);
        GeneralCollisionData<2,25> prams(cfg, tgv, scaling, viscosity, d2q25h,
                                        cs2 , dt);

        prams.velocity[0]=velocities[0];
        prams.velocity[1]=velocities[1];
        prams.density = rho;
        prams.temperature = 1.1;
        prams.H3 = calculateH3<2,25>(prams);
        prams.H4 = calculateH4<2,25>(prams);

        eq.calc(feq,prams);

        /*for (int i =0;i<9;i++)
        {
            cout << feq[i] << endl;
        }*/

        //cout << "Dichte v:" << prams.density << endl;
        //cout << "Geschwindigkeit vorher:" << prams.velocity[0] << endl;
        //cout << "Geschwindigkeit vorher:" << prams.velocity[1] << endl;
        std::array<double,2> v_post;
        calculateVelocity<2,25>(feq,v_post,calculateDensity<25>(feq),prams);

        BOOST_CHECK_CLOSE(v_post[0],prams.velocity[0],10e-6);
        BOOST_CHECK_CLOSE(v_post[1],prams.velocity[1],10e-6);





    } // test_case

    BOOST_AUTO_TEST_CASE(Equilibrium_D3Q77_test) {



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

        double velocities[3]={0.1,0.1,0.1};

        QuarticEquilibrium<3,27> eq;
        std::array<double,27> feq;
        D3Q27 d3q27(1.0);
        double cs2=1./3.;
        double scaling = 1.0;
        double dt = 0.1;

        double viscosity = 1.0;


        SolverConfiguration cfg;
        TaylorGreenVortex3D tgv(0.1,1);
        GeneralCollisionData<3,27> prams(cfg, tgv, scaling, viscosity, d3q27,
                                         cs2 , dt);

        prams.velocity[0]=velocities[0];
        prams.velocity[1]=velocities[1];
        prams.velocity[2]=velocities[2];
        prams.density = rho;
        prams.temperature = 1.3;
        prams.H3 = calculateH3<3,27>(prams);
        prams.H4 = calculateH4<3,27>(prams);

        eq.calc(feq,prams);

        /*for (int i =0;i<9;i++)
        {
            cout << feq[i] << endl;
        }*/

        //cout << "Dichte v:" << prams.density << endl;
        //cout << "Geschwindigkeit vorher:" << prams.velocity[0] << endl;
        //cout << "Geschwindigkeit vorher:" << prams.velocity[1] << endl;
        for (auto& i : feq)
        {
            cout << i << endl;
        }
        std::array<double,3> v_post;
        double rho_post = calculateDensity<27>(feq);
        calculateVelocity<3,27>(feq,v_post,calculateDensity<27>(feq),prams);

        BOOST_CHECK_CLOSE(rho_post,rho,10e-6);

        BOOST_CHECK_CLOSE(v_post[0],prams.velocity[0],10e-6);
        BOOST_CHECK_CLOSE(v_post[1],prams.velocity[1],10e-6);
        BOOST_CHECK_CLOSE(v_post[2],prams.velocity[2],10e-6);






    } // test_case


BOOST_AUTO_TEST_SUITE_END()




