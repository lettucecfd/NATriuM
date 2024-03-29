/**
 * @file VelocityNeqBounceBack_test.cpp
 * @short
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */


#include "natrium/boundaries/VelocityNeqBounceBack.h"

#include "./WallFixture.h"
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


using namespace natrium;

BOOST_AUTO_TEST_SUITE(VelocityNeqBounceBack2D_DG_test)



BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack2D_DG_Construction_test) {
	pout << "VelocityNeqBounceBack_DG_Construction_test..." << endl;

	BOOST_CHECK_NO_THROW(
			VelocityNeqBounceBack<2> mlBound1(0, boost::make_shared<BoundaryTestVelocity>()));
	numeric_vector U(2);
	BOOST_CHECK_NO_THROW(VelocityNeqBounceBack<2> mlBound2(0, U); );

	pout << "done" << endl;
} /*VelocityNeqBounceBack2D_DG_Construction_test */

BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack2D_DG_SparsityPattern_test) {
	pout << "VelocityNeqBounceBack_DG_SparsityPattern_test..." << endl;

	// The incoming particle distributions at the boundary must be affected by the opposite outgoing ones
	// This means that diagonal entries must exist for the boundary dofs
	// for all blocks (I,J) (I for incoming and J for their opposites)
	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
			WallTestDomain2D>(1);
	SEDGMinLee<2> advector(*problem, 2,
			boost::make_shared<D2Q9>());
	advector.setupDoFs();
	advector.reassemble();
	vector<bool> isBoundary(advector.getNumberOfDoFs());
	for (size_t i = 0; i < advector.getNumberOfDoFs(); i++) {
		std::set<dealii::types::boundary_id> boundaryIndicators;
		// left boundary
		boundaryIndicators.insert(0);
		/*		dealii::DoFTools::extract_dofs_with_support_on_boundary(*(advector.getDoFHandler()),
		 dealii::ComponentMask(), isBoundary, boundaryIndicators);
		 if (isBoundary.at(i)) {
		 //pout << i << endl;
		 // note that block 0 refers to f_1 and so on
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(4, 6).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(0, 2).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(5, 7).exists(i,
		 i));
		 }
		 // right boundary
		 boundaryIndicators.clear();
		 boundaryIndicators.insert(1);
		 dealii::DoFTools::extract_dofs_with_support_on_boundary(*(advector.getDoFHandler()),
		 dealii::ComponentMask(), isBoundary, boundaryIndicators);
		 if (isBoundary.at(i)) {
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(6, 4).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(2, 0).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(7, 5).exists(i,
		 i));
		 }
		 // bottom boundary
		 boundaryIndicators.clear();
		 boundaryIndicators.insert(2);
		 dealii::DoFTools::extract_dofs_with_support_on_boundary(*(advector.getDoFHandler()),
		 dealii::ComponentMask(), isBoundary, boundaryIndicators);
		 if (isBoundary.at(i)) {
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(4, 6).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(1, 3).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(5, 7).exists(i,
		 i));
		 }
		 // top boundary
		 boundaryIndicators.clear();
		 boundaryIndicators.insert(3);
		 dealii::DoFTools::extract_dofs_with_support_on_boundary(*(advector.getDoFHandler()),
		 dealii::ComponentMask(), isBoundary, boundaryIndicators);
		 if (isBoundary.at(i)) {
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(6, 4).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(3, 1).exists(i,
		 i));
		 BOOST_CHECK(
		 advector.getBlockSparsityPattern().block(7, 5).exists(i,
		 i));
		 }
		 */
	}

	pout << "done" << endl;
} /* VelocityNeqBounceBack2D_DG_SparsityPattern_test */

BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack_DG_MassConservation_test) {
	pout << "VelocityNeqBounceBack_DG_MassConservation_test..." << endl;

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

	CFDSolver<2> solver(configuration, problem);
	solver.run();

	// check mass conservation
	double mass = solver.getDensity().l1_norm();
	mass /= solver.getNumberOfDoFs();
	BOOST_CHECK_SMALL(mass - 1.0, 1e-10);

	pout << "done" << endl;
} /*VelocityNeqBounceBack2D_DG_MassConservation_test */


BOOST_AUTO_TEST_CASE(VelocityNeqBounceBack2D_DG_BoundaryVelocity_test) {
	pout << "VelocityNeqBounceBack2D_DG_BoundaryVelocity_test..." << endl;

	WallTest w(SEDG);

	// Check velocity in domain
	BOOST_CHECK_SMALL(w.error_u, 1e-3);
	BOOST_CHECK_SMALL(w.error_v, 1e-5);


	pout << "done" << endl;
} /* VelocityNeqBounceBack2D_DG_BoundaryVelocity_test */




BOOST_AUTO_TEST_SUITE_END() /*VelocityNeqBounceBack2D_test*/

