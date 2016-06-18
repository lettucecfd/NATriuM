/**
 * @file LinearBoundaryRhoU2D_test.cpp
 * @short
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <natrium/problemdescription/LinearBoundaryRhoU.h>
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
#include "natrium/advection/SemiLagrangianBoundaryDoFHandler.h"
#include "natrium/advection/SemiLagrangianVectorReferenceTypes.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "WallTestDomain2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(LinearBoundaryRhoU2D_test)

class BoundaryTestDensity: public dealii::Function<2> {
public:
	virtual double value(const dealii::Point<2> &, const unsigned int) const {
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

BOOST_AUTO_TEST_CASE(LinearBoundaryRhoU2D_Construction_test) {
	pout << "LinearBoundaryRhoU2D_Construction_test..." << endl;

	BOOST_CHECK_NO_THROW(
			LinearBoundaryRhoU<2> mlBound1(0, boost::make_shared<BoundaryTestDensity>(), boost::make_shared<BoundaryTestVelocity>()));
	numeric_vector U(2);
	BOOST_CHECK_NO_THROW(LinearBoundaryRhoU<2> mlBound2(0, U); );

	pout << "done" << endl;
} /*LinearBoundaryRhoU2D_Construction_test */

BOOST_AUTO_TEST_CASE(LinearBoundaryRhoU2D_SparsityPattern_test) {
	pout << "LinearBoundaryRhoU2D_SparsityPattern_test..." << endl;

	// The incoming particle distributions at the boundary must be affected by the opposite outgoing ones
	// This means that diagonal entries must exist for the boundary dofs
	// for all blocks (I,J) (I for incoming and J for their opposites)
	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
			WallTestDomain2D>(1);
	SEDGMinLee<2> advector(problem->getMesh(), problem->getBoundaries(), 2,
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
} /* LinearBoundaryRhoU2D_SparsityPattern_test */

BOOST_AUTO_TEST_CASE(LinearBoundaryRhoU2D_MassConservation_test) {
	pout << "LinearBoundaryRhoU2D_MassConservation_test..." << endl;

	// make problem and solver
	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
			WallTestDomain2D>(1);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("/tmp");
	configuration->setSwitchOutputOff(true);
	configuration->setNumberOfTimeSteps(100);
	configuration->setSedgOrderOfFiniteElement(1);
	configuration->setCFL(0.4);

	CFDSolver<2> solver(configuration, problem);

	solver.run();

	// check mass conservation
	double mass = solver.getDensity().l1_norm();
	mass /= solver.getNumberOfDoFs();
	BOOST_CHECK_SMALL(mass - 1.0, 1e-10);

	pout << "done" << endl;
} /*LinearBoundaryRhoU2D_MassConservation_test */

BOOST_AUTO_TEST_CASE(LinearBoundaryRhoU2D_BoundaryVelocity_test) {
	pout << "LinearBoundaryRhoU2D_BoundaryVelocity_test..." << endl;

	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
			WallTestDomain2D>(1);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("/tmp");
	configuration->setSwitchOutputOff(true);
	configuration->setNumberOfTimeSteps(1);
	configuration->setSedgOrderOfFiniteElement(9);
	configuration->setCFL(0.4);

	CFDSolver<2> solver(configuration, problem);
	solver.run();

	// Check boundary velocity
	std::set<dealii::types::boundary_id> boundaryIndicators;
	boundaryIndicators.insert(3);
	std::vector<bool> isBoundary(solver.getNumberOfDoFs());
	DealIIExtensions::extract_dofs_with_support_on_boundary(
			*(solver.getAdvectionOperator()->getDoFHandler()),
			dealii::ComponentMask(), isBoundary, boundaryIndicators);
	/*
	 // TODO good test
	 for (size_t i = 0; i < isBoundary.size(); i++) {
	 if (isBoundary.at(i)) {
	 pout << i << endl;
	 BOOST_CHECK_CLOSE(
	 solver.getVelocity().at(0)(i),
	 0.01, 0.00001);
	 BOOST_CHECK_SMALL(solver.getVelocity().at(1)(i), 1e-7);
	 }
	 }
	 */

	pout << "done" << endl;
} /* LinearBoundaryRhoU2D_BoundaryVelocity_test */

BOOST_AUTO_TEST_CASE(LinearBoundaryRhoU2D_makeIncomingDirections_test) {
	pout << "LinearBoundaryRhoU2D_makeIncomingDirections_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();

	dealii::DoFHandler<2> dof_handler(*periodic.getMesh());
	dealii::FE_DGQArbitraryNodes<2> fe(
			dealii::QGaussLobatto<1>(fe_order + 1));
	dof_handler.distribute_dofs(fe);

	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			dof_handler.begin_active();
	BoundaryHit<2> hit(dealii::Point<2>(0.0, 0.0), 0.0, dealii::Tensor<1, 2>(),
			*(periodic.getBoundaries()->getBoundary(0)), cell, GeneralizedDoF(false,0,1));
	D2Q9 d2q9;
	numeric_vector ub1(2);
	LinearBoundaryRhoU<2> boundary1(0, ub1);
	boundary1.makeIncomingDirections(hit, d2q9);

	BOOST_CHECK_EQUAL(hit.in.size(), size_t(1));
	BOOST_CHECK_EQUAL(hit.in.at(0).getAlpha(), size_t(3));

	dof_handler.clear();
	pout << "done" << endl;
} /* LinearBoundaryRhoU2D_makeIncomingDirections_test */


BOOST_AUTO_TEST_CASE(LinearBoundaryRhoU2D_calculate_test) {
	pout << "LinearBoundaryRhoU2D_calculate_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();

	dealii::DoFHandler<2> dof_handler(*periodic.getMesh());
	dealii::FE_DGQArbitraryNodes<2> fe(
			dealii::QGaussLobatto<1>(fe_order + 1));
	dof_handler.distribute_dofs(fe);
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			dof_handler.begin_active();
	vector<distributed_vector> v1;
	for (size_t i = 0; i < 9; i++){
		distributed_vector tmp(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
		v1.push_back(tmp);
	}
	DistributionFunctions fold (v1);

	vector<distributed_vector> v2;
	for (size_t i = 0; i < 9; i++){
		distributed_vector tmp(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
		v2.push_back(tmp);
	}
	DistributionFunctions fnew (v2);

	SecondaryBoundaryDoFVector sbdv(dof_handler.locally_owned_dofs());

	SemiLagrangianVectorAccess generalized_f(fold, fnew, sbdv);

	GeneralizedDoF dof_out(false, 0,1);
	BoundaryHit<2> hit(dealii::Point<2>(0.0, 0.0), 0.0, dealii::Tensor<1, 2>(),
			*(periodic.getBoundaries()->getBoundary(0)), cell, dof_out);

	// simple BB
	D2Q9 d2q9;
	numeric_vector ub1(2);
	LinearBoundaryRhoU<2> boundary1(0, ub1);
	boundary1.makeIncomingDirections(hit, d2q9);
	GeneralizedDoF dof_in(false, 0, 3);
	hit.in.push_back(dof_in);
	boundary1.calculate(hit, d2q9, 0, generalized_f);
	BOOST_CHECK_SMALL(generalized_f[dof_out] - 1.0, 1e-15);


	// velocity BB
	generalized_f[dof_out] = 1.0;
	BoundaryHit<2> hit2(dealii::Point<2>(0.0, 0.0), 0.0, dealii::Tensor<1, 2>(),
			*(periodic.getBoundaries()->getBoundary(0)), cell, dof_out);
	numeric_vector ub2(2);
	ub2(0) = 2.5;
	LinearBoundaryRhoU<2> boundary2(0, ub2);
	boundary2.makeIncomingDirections(hit2, d2q9);
	hit2.in.push_back(dof_in);
	boundary2.calculate(hit2, d2q9, 0, generalized_f);
	double expected = 1.0 + 2.0 * 1.0 / 9.0 * 2.5 / (1.0/3.0);
	// f_opposite + 2 * stencil.getWeight(boundary_hit.outgoingDirection) * 1
	// * (ea * velocity) / stencil.getSpeedOfSoundSquare()

	BOOST_CHECK_SMALL(generalized_f[dof_out] - expected, 1e-15);

	dof_handler.clear();
	pout << "done" << endl;
} /* LinearBoundaryRhoU2D_calculate_test */

BOOST_AUTO_TEST_SUITE_END() /*LinearBoundaryRhoU2D_test*/

} /* namespace natrium */
