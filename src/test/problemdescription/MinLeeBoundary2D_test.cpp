/**
 * @file MinLeeBoundary2D_test.cpp
 * @short
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/problemdescription/MinLeeBoundary.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/base/function.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/fe/component_mask.h"
#include "deal.II/base/index_set.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/advection/SEDGMinLee.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "WallTestDomain2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(MinLeeBoundary2D_test)

class BoundaryTestDensity: public dealii::Function<2> {
public:
	virtual double value(const dealii::Point<2> &,
			const unsigned int ) const {
		return 1;
	}
};
class BoundaryTestVelocity: public dealii::Function<2> {
public:
	virtual void vector_value(const dealii::Point<2> &,
			dealii::Vector<double> &values) const {
		values(0) = 0;
		values(1) = 0;
	}
};

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_Construction_test) {
	pout << "MinLeeBoundary2D_Construction_test..." << endl;

	BOOST_CHECK_NO_THROW(
			MinLeeBoundary<2> mlBound1(0, make_shared<BoundaryTestDensity>(), make_shared<BoundaryTestVelocity>()));
	numeric_vector U(2);
	BOOST_CHECK_NO_THROW(MinLeeBoundary<2> mlBound2(0, U); );

	pout << "done" << endl;
} /*MinLeeBoundary2D_Construction_test */

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_SparsityPattern_test) {
	pout << "MinLeeBoundary2D_SparsityPattern_test..." << endl;

	// The incoming particle distributions at the boundary must be affected by the opposite outgoing ones
	// This means that diagonal entries must exist for the boundary dofs
	// for all blocks (I,J) (I for incoming and J for their opposites)
	shared_ptr<ProblemDescription<2> > problem = make_shared<WallTestDomain2D>(
			1);
	SEDGMinLee<2> advector(problem->getMesh(), problem->getBoundaries(), 2,
			make_shared<D2Q9>());
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
} /* MinLeeBoundary2D_SparsityPattern_test */

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_MassConservation_test) {
	pout << "MinLeeBoundary2D_MassConservation_test..." << endl;

	// make problem and solver
	shared_ptr<ProblemDescription<2> > problem = make_shared<WallTestDomain2D>(
			1);
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("/tmp");
	configuration->setSwitchOutputOff(true);
	configuration->setNumberOfTimeSteps(100);
	configuration->setSedgOrderOfFiniteElement(1);
	configuration->setTimeStepSize(0.01);

	CFDSolver<2> solver(configuration, problem);

	solver.run();

	// check mass conservation
	double mass = 0.0;
	for (size_t i = 0; i < solver.getNumberOfDoFs(); i++) {
		mass += solver.getDensity()(i);
	}
	mass /= solver.getNumberOfDoFs();
	BOOST_CHECK_SMALL(mass - 1.0, 1e-10);

	pout << "done" << endl;
} /*MinLeeBoundary2D_MassConservation_test */

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_BoundaryVelocity_test) {
	pout << "MinLeeBoundary2D_BoundaryVelocity_test..." << endl;

	shared_ptr<ProblemDescription<2> > problem = make_shared<WallTestDomain2D>(
			1);
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("/tmp");
	configuration->setSwitchOutputOff(true);
	configuration->setNumberOfTimeSteps(1);
	configuration->setSedgOrderOfFiniteElement(9);
	configuration->setTimeStepSize(0.0001);

	CFDSolver<2> solver(configuration, problem);
	solver.run();

	// Check boundary velocity
	std::set<dealii::types::boundary_id> boundaryIndicators;
	boundaryIndicators.insert(3);
	dealii::IndexSet isBoundary;
	dealii::DoFTools::extract_boundary_dofs(
			*(solver.getAdvectionOperator()->getDoFHandler()),
			dealii::ComponentMask(), isBoundary);	//, boundaryIndicators);

	for (size_t i = 0; i < isBoundary.n_elements(); i++) {
		pout << i << endl;
		BOOST_CHECK_CLOSE(
				solver.getVelocity().at(0)(isBoundary.nth_index_in_set(i)),
				0.01, 0.00001);
		BOOST_CHECK_SMALL(solver.getVelocity().at(1)(i), 1e-7);

	}

	pout << "done" << endl;
} /* MinLeeBoundary2D_BoundaryVelocity_test */

BOOST_AUTO_TEST_SUITE_END() /*MinLeeBoundary2D_test*/

} /* namespace natrium */
