

/**
 * @file VelocityNeqBounceBack_SL_test.cpp
 * @short
 * @date 15.02.2018
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */


#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include "boost/test/unit_test.hpp"

#include "deal.II/base/function.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/fe/component_mask.h"
#include "deal.II/base/index_set.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/advection/SEDGMinLee.h"
#include "natrium/advection/SemiLagrangian.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "../problemdescription/WallTestDomain2D.h"
#include "WallFixture.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(VelocityNeqBounceBack2D_SL_test)


//BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack2D_SL_Construction_test) {
//	pout << "VelocityNeqBounceBack2D_SL_Construction_test..." << endl;
//
//	dealii::Tensor<1,2> U;
//	BOOST_CHECK_NO_THROW(
//			VelocityNeqBounceBack<2> mlBound1(0,  U));
//	/*numeric_vector U(2);
//	BOOST_CHECK_NO_THROW(VelocityNeqBounceBack<2> mlBound2(0, U); );*/
//
//	pout << "done" << endl;
//} /*VelocityNeqBounceBack2D_SL_Construction_test */



BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack2D_SL_MassConservation_test) {
	pout << "VelocityNeqBounceBack2D_SL_MassConservation_test..." << endl;

	// make problem and solver
	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
			WallTestDomain2D>(1);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("/tmp");
	configuration->setSwitchOutputOff(true);
	configuration->setNumberOfTimeSteps(10);
	configuration->setSedgOrderOfFiniteElement(1);
	configuration->setCFL(0.4);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);

	CFDSolver<2> solver(configuration, problem);
	solver.run();

	// check mass conservation
	double mass = solver.getDensity().l1_norm();
	mass /= solver.getNumberOfDoFs();
	BOOST_CHECK_SMALL(mass - 1.0, 1e-10);

	pout << "done" << endl;
} /*VelocityNeqBounceBack2D_SL_MassConservation_test */


BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack2D_SL_BoundaryVelocity_test) {
	pout << "VelocityNeqBounceBack2D_SL_BoundaryVelocity_test..." << endl;

	WallTest w(SEMI_LAGRANGIAN);

	// Check velocity in domain
	BOOST_CHECK_SMALL(w.error_u, 1e-3);
	BOOST_CHECK_SMALL(w.error_v, 1e-5);


	pout << "done" << endl;
} /* VelocityNeqBounceBack2D_SL_BoundaryVelocity_test */





BOOST_AUTO_TEST_SUITE_END() /*VelocityNeqBounceBack2D_SL_test*/


