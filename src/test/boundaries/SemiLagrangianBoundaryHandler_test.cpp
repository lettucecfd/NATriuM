/*
 * SemiLagrangianBoundaryHandler_test.cpp
 *
 *  Created on: 17.11.2016
 *      Author: akraem3m
 */

#include "natrium/boundaries/SemiLagrangianBoundaryHandler.h"

#include "boost/test/unit_test.hpp"

#include "natrium/problemdescription/BoundaryCollection.h"
#include "natrium/advection/SemiLagrangian.h"
#include "natrium/advection/SemiLagrangianTools.h"
#include "natrium/benchmarks/CouetteFlow2D.h"
#include "natrium/benchmarks/CouetteFlow3D.h"
//#include "natrium/benchmarks/CouetteFlowGrad2D.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/boundaries/SLFirstOrderBounceBack.h"
#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/benchmarks/PoiseuilleFlow2D.h"

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(SemiLagrangianBoundaryHandler_test)

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_Construction_test) {
	pout << "SemiLagrangianBoundaryHandler_Construction_test..." << endl;

	double viscosity = 1e-1;
	double u0 = 0.1;
	size_t ref_level = 3;
	size_t fe_order = 2;

	// D2Q9
	boost::shared_ptr<ProblemDescription<2> > couette = boost::make_shared<
			CouetteFlow2D>(viscosity, u0, ref_level);
	couette->refineAndTransform();
	SemiLagrangian<2> streaming(*couette, fe_order, boost::make_shared<D2Q9>(),
			0.001);
	SemiLagrangianBoundaryHandler<2> bh(streaming);

	// D3Q19
	boost::shared_ptr<ProblemDescription<3> > couette3d = boost::make_shared<
			CouetteFlow3D>(viscosity, u0, ref_level);
	couette->refineAndTransform();
	SemiLagrangian<3> streaming_d3q19(*couette3d, fe_order,
			boost::make_shared<D3Q19>(), 0.001);
	SemiLagrangianBoundaryHandler<3> bh_d3q19(streaming_d3q19);

	// D3Q15
	SemiLagrangian<3> streaming_d3q15(*couette3d, fe_order,
			boost::make_shared<D3Q15>(), 0.001);
	SemiLagrangianBoundaryHandler<3> bh_d3q15(streaming_d3q15);

	// D3Q27
	SemiLagrangian<3> streaming_d3q27(*couette3d, fe_order,
			boost::make_shared<D3Q27>(), 0.001);
	SemiLagrangianBoundaryHandler<3> bh_d3q27(streaming_d3q27);

	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_Construction_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_addHit_test) {
	pout << "SemiLagrangianBoundaryHandler_addHit_test..." << endl;

	D2Q9 d2q9;

	double viscosity = 1e-1;
	double u0 = 0.1;
	const double dt = 1e-4;
	size_t ref_level = 3;
	size_t fe_order = 2;
	CouetteFlow2D couette(viscosity, u0, ref_level);
	couette.refineAndTransform();
	SemiLagrangian<2> streaming(couette, fe_order, boost::make_shared<D2Q9>(),
			dt);
	streaming.setupDoFs();

	BoundaryCollection<2> bc2;
	dealii::Point<2> departure(-1, -1);
	dealii::Point<2> p(0, 0);
	LagrangianPathTracker<2> tracker(1, 1, 1, departure, p,
			streaming.getDoFHandler()->begin_active());
	assert(streaming.getStencil());
	SemiLagrangianBoundaryHandler<2> bh(streaming);
	bh.addHit(tracker, 2);
	BOOST_CHECK(1 == bh.n_cells());

	bh.addHit(tracker, 2);
	BOOST_CHECK(1 == bh.n_cells());

	LagrangianPathTracker<2> tracker2(1, 1, 1, departure, p,
			++(streaming.getDoFHandler()->begin_active()));
	bh.addHit(tracker2, 2);
	BOOST_CHECK(2 == bh.n_cells());

	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_addHit_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_DataStructures_test) {
	pout << "SemiLagrangianBoundaryHandler_DataStructures_test..." << endl;

	double viscosity = 1e-1;
	double u0 = 0.1;
	size_t ref_level = 3;
	size_t fe_order = 2;
	boost::shared_ptr<ProblemDescription<2> > couette = boost::make_shared<
			CouetteFlow2D>(viscosity, u0, ref_level);
	couette->refineAndTransform();
	SemiLagrangian<2> streaming(*couette, fe_order, boost::make_shared<D2Q9>(),
			0.001);
	streaming.setupDoFs();
	streaming.reassemble();

	// small time step => only support hits
	// precondition: all processes have an equal number of cells and boundary cells
	const size_t expected_n_cells = (pow(2, ref_level) * 2)
			/ dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	BOOST_CHECK_EQUAL(expected_n_cells,
			streaming.getBoundaryHandler().n_cells());
	const size_t n_undefined_directions = 3;
	const size_t n_points = pow(2, ref_level) * 2 * (fe_order) + 2;
	const size_t expected_n_hits = n_points * n_undefined_directions
			/ dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	BOOST_CHECK_EQUAL(expected_n_hits, streaming.getBoundaryHandler().n_hits());

	// TODO: finalize
	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_DataStructures_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_offSupportHits_test) {
	pout << "SemiLagrangianBoundaryHandler_offSupportHits_test..." << endl;

	double viscosity = 1e-1;
	double u0 = 0.1;
	size_t ref_level = 3;
	size_t fe_order = 2;
	boost::shared_ptr<ProblemDescription<2> > couette = boost::make_shared<
			CouetteFlow2D>(viscosity, u0, ref_level);
	couette->refineAndTransform();
	SemiLagrangian<2> streaming(*couette, fe_order, boost::make_shared<D2Q9>(),
			0.063);
	// time step chosen > 0.0625, i.e., second layer of lattice nodes belong to boundary
	streaming.setupDoFs();
	streaming.reassemble();

	// small time step => only support hits
	// precondition: all processes have an equal number of cells and boundary cells
	const size_t expected_n_cells = (pow(2, ref_level) * 2)
			/ dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	BOOST_CHECK_EQUAL(expected_n_cells,
			streaming.getBoundaryHandler().n_cells());
	const size_t n_undefined_directions = 3;
	// twice as many points as in the previous test:
	const size_t n_points = 2 * (pow(2, ref_level) * 2 * (fe_order) + 2);
	const size_t expected_n_hits = n_points * n_undefined_directions
			/ dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	BOOST_CHECK_EQUAL(expected_n_hits, streaming.getBoundaryHandler().n_hits());

	// TODO: finalize
	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_offSupportHits_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_PoiseuilleBB_test) {
	pout << "SemiLagrangianBoundaryHandler_PoiseuilleBB_test..." << endl;

	// setup test case
	const double CFL = 1.5;
	const double Re = 10;
	const double u_bulk = 0.0001 / 1.5; //1.0;
	const double height = 1.0;
	const double length = 2.0;
	const double orderOfFiniteElement = 2;
	const double Ma = 0.1;
	const double refinement_level = 2;
	bool is_periodic = true;

	/// create CFD problem
	double viscosity = u_bulk * height / Re;
	const double scaling = sqrt(3) * 1.5 * u_bulk / Ma;

	/// setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setUserInteraction(false);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(CFL);
	configuration->setConvergenceThreshold(1e-6);
	configuration->setForcingScheme(SHIFTING_VELOCITY);

	// make solver object and run simulation
	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
			PoiseuilleFlow2D>(viscosity, refinement_level, u_bulk, height,
			length, is_periodic);
	// overwrite boundary conditions
	boost::shared_ptr<BoundaryCollection<2> > bc = boost::make_shared<
			BoundaryCollection<2> >();
	// we have to take care that we do not create a periodic boundary twice
	bc->addBoundary(problem->getBoundaries()->getPeriodicBoundary(0));
	bc->addBoundary(boost::make_shared<SLFirstOrderBounceBack<2> >(2));
	bc->addBoundary(boost::make_shared<SLFirstOrderBounceBack<2> >(3));
	problem->setBoundaries(bc);

	CFDSolver<2> solver(configuration, problem);

	//BOOST_CHECK(solver.getVelocity().at(0).norm_sqr() < 1e-10);
	solver.run();
	solver.getSolverStats()->update();
	BOOST_CHECK_GT(solver.getSolverStats()->getMeanVelocityX(), 0.9*u_bulk);
	BOOST_CHECK_LT(solver.getSolverStats()->getMeanVelocityX(), 1.1*u_bulk);

	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_PoiseuilleBB_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_Couette_test) {
	pout << "SemiLagrangianBoundaryHandler_Couette_test..." << endl;

	double viscosity = 1e-1;
	double u0 = 0.1;
	size_t ref_level = 3;
	boost::shared_ptr<ProblemDescription<2> > couette = boost::make_shared<
			CouetteFlow2D>(viscosity, u0, ref_level);

	// TODO: finalize
	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_Couette_test */

BOOST_AUTO_TEST_SUITE_END()
