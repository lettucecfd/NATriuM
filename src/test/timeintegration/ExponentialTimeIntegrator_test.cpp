/**
 * @file ExponentialTimeIntegrator_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <iostream>

#include "boost/test/unit_test.hpp"
#include "deal.II/lac/sparsity_pattern.h"
#include "deal.II/lac/compressed_sparsity_pattern.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/timeintegration/ExponentialTimeIntegrator.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(ExponentialTimeIntegrator_test)

BOOST_AUTO_TEST_CASE(ExponentialTimeIntegrator_test) {


	cout << "ExponentialTimeIntegrator_test..." << endl;

	// solve ODE F(f) = lambda*f
	// with initial value f0 = 1
	// analytic solution: f(x) = e^(lamda*x)
	const double dt = 0.01;
	const size_t numberOfSteps = 100;
	distributed_vector f(10);

	f(0)=-3;
	f(1)=0;
	f(2)=0;
	f(3)=4;
	f(4)=-4;
	f(5)=-2;
	f(6)=-3;
	f(7)=-2;
	f(8)=-5;
	f(9)=-5;
	distributed_vector g(10);

	g(0)=2;
	g(1)=0;
	g(2)=-4;
	g(3)=4;
	g(4)=0;
	g(5)=1;
	g(6)=-4;
	g(7)=-1;
	g(8)=-5;
	g(9)=-1;


	// build the 1x1 matrix [[lambda]]
	dealii::CompressedSparsityPattern compressedSparsityPattern(10,10);
	for (int i = 0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			compressedSparsityPattern.add(i,j);
		}
	}
	dealii::SparsityPattern sparsityPattern;
	sparsityPattern.copy_from(compressedSparsityPattern);
	distributed_sparse_matrix A;
	A.reinit(sparsityPattern);
	A.set(0,0,-1);
	A.set(0,1,1);
	A.set(0,2,-4);
	A.set(0,3,-2);
	A.set(0,4,-4);
	A.set(0,5,-1);
	A.set(0,6,4);
	A.set(0,7,-3);
	A.set(0,8,-4);
	A.set(0,9,5);
	A.set(1,0,-3);
	A.set(1,1,-3);
	A.set(1,2,0);
	A.set(1,3,1);
	A.set(1,4,5);
	A.set(1,5,2);
	A.set(1,6,1);
	A.set(1,7,1);
	A.set(1,8,-4);
	A.set(1,9,4);
	A.set(2,0,0);
	A.set(2,1,-3);
	A.set(2,2,-3);
	A.set(2,3,2);
	A.set(2,4,2);
	A.set(2,5,0);
	A.set(2,6,-2);
	A.set(2,7,5);
	A.set(2,8,0);
	A.set(2,9,-4);
	A.set(3,0,1);
	A.set(3,1,-2);
	A.set(3,2,4);
	A.set(3,3,2);
	A.set(3,4,-3);
	A.set(3,5,-2);
	A.set(3,6,1);
	A.set(3,7,4);
	A.set(3,8,-1);
	A.set(3,9,0);
	A.set(4,0,4);
	A.set(4,1,-5);
	A.set(4,2,2);
	A.set(4,3,-4);
	A.set(4,4,2);
	A.set(4,5,-4);
	A.set(4,6,-4);
	A.set(4,7,2);
	A.set(4,8,-2);
	A.set(4,9,1);
	A.set(5,0,5);
	A.set(5,1,-4);
	A.set(5,2,-2);
	A.set(5,3,0);
	A.set(5,4,-1);
	A.set(5,5,4);
	A.set(5,6,2);
	A.set(5,7,-4);
	A.set(5,8,5);
	A.set(5,9,1);
	A.set(6,0,-3);
	A.set(6,1,4);
	A.set(6,2,4);
	A.set(6,3,1);
	A.set(6,4,0);
	A.set(6,5,4);
	A.set(6,6,4);
	A.set(6,7,-1);
	A.set(6,8,-1);
	A.set(6,9,3);
	A.set(7,0,-3);
	A.set(7,1,-1);
	A.set(7,2,-1);
	A.set(7,3,0);
	A.set(7,4,5);
	A.set(7,5,-1);
	A.set(7,6,-2);
	A.set(7,7,2);
	A.set(7,8,5);
	A.set(7,9,4);
	A.set(8,0,2);
	A.set(8,1,0);
	A.set(8,2,1);
	A.set(8,3,-4);
	A.set(8,4,-1);
	A.set(8,5,-3);
	A.set(8,6,3);
	A.set(8,7,2);
	A.set(8,8,1);
	A.set(8,9,1);
	A.set(9,0,4);
	A.set(9,1,4);
	A.set(9,2,-3);
	A.set(9,3,3);
	A.set(9,4,-1);
	A.set(9,5,-5);
	A.set(9,6,-4);
	A.set(9,7,1);
	A.set(9,8,0);
	A.set(9,9,2);
	distributed_vector z(10);

	z(0)=-2.251746e+03;
	z(1)=3.070874e+03;
	z(2)=2.294584e+03;
	z(3)=6.848658e+03;
	z(4)=1.282751e+03;
	z(5)=-8.276442e+03;
	z(6)=2.004676e+03;
	z(7)=1.151836e+04;
	z(8)=4.635269e+03;
	z(9)=8.616738e+03;

	// initialize Runge-Kutta 5
	ExponentialTimeIntegrator<distributed_sparse_matrix, distributed_vector> ETI(dt, 1); //

	double t = 0;
	for (size_t i = 0; i < numberOfSteps; i++){
		t += dt;
		ETI.step(f,A,g);



	}

	double error = fabs(f(5))-fabs(z(5));

	BOOST_ASSERT(error < 1e-2);

	cout << "done." << endl;


}
	BOOST_AUTO_TEST_SUITE_END()





} /* namespace natrium */
