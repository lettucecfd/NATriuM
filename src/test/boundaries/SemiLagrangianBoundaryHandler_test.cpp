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

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(SemiLagrangianBoundaryHandler_test)

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_Construction_test) {
	pout << "SemiLagrangianBoundaryHandler_Construction_test..." << endl;

	D2Q9 d2q9;
	D3Q19 d3q19;
	D3Q15 d3q15;
	D3Q27 d3q27;
	BoundaryCollection<2> bc2;
	BoundaryCollection<3> bc3;
	BOOST_CHECK_NO_THROW(SemiLagrangianBoundaryHandler<2>(1e-4, d2q9, bc2));
	BOOST_CHECK_NO_THROW(SemiLagrangianBoundaryHandler<3>(1e-4, d3q19, bc3));
	BOOST_CHECK_NO_THROW(SemiLagrangianBoundaryHandler<3>(1e-4, d3q15, bc3));
	BOOST_CHECK_NO_THROW(SemiLagrangianBoundaryHandler<3>(1e-4, d3q27, bc3));

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
	CouetteFlowGrad2D couette(viscosity, u0, ref_level);
	couette.refineAndTransform();
	SemiLagrangian<2> streaming(couette.getMesh(), couette.getBoundaries(),
			fe_order, boost::make_shared<D2Q9>(), dt);
	streaming.setupDoFs();

	BoundaryCollection<2> bc2;
	dealii::Point<2> departure(-1, -1);
	dealii::Point<2> p(0,0);
	LagrangianPathTracker<2> tracker(1, 1, 1, departure, p,
			streaming.getDoFHandler()->begin_active());
	assert (streaming.getStencil());
	SemiLagrangianBoundaryHandler<2> bh(dt, *streaming.getStencil(), *couette.getBoundaries());
	bh.addHit(tracker, 2, streaming);
	BOOST_CHECK(1 == bh.n_cells());

	bh.addHit(tracker, 2, streaming);
	BOOST_CHECK(1 == bh.n_cells());

	LagrangianPathTracker<2> tracker2(1, 1, 1, departure, p,
			++(streaming.getDoFHandler()->begin_active()));
	bh.addHit(tracker2, 2, streaming);
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
			CouetteFlowGrad2D>(viscosity, u0, ref_level);
	couette->refineAndTransform();
	SemiLagrangian<2> streaming(couette->getMesh(), couette->getBoundaries(),
			fe_order, boost::make_shared<D2Q9>(), 0.001);
	streaming.setupDoFs();
	streaming.reassemble();

	// small time step => only support hits
	// precondition: all processes have an equal number of cells and boundary cells
	const size_t expected_n_cells = (pow(2, ref_level) * 2)
			/ dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	BOOST_CHECK_EQUAL(expected_n_cells,
			streaming.getBoundaryHandler().n_cells());
	const size_t n_undefined_directions = 3;
	const size_t expected_n_hits = (pow(2, ref_level) * 2 * (fe_order + 1) * n_undefined_directions)
			/ dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); // minus corners
	BOOST_CHECK_EQUAL(expected_n_hits, streaming.getBoundaryHandler().n_hits());
	BOOST_CHECK_EQUAL(expected_n_hits,
			streaming.getBoundaryHandler().n_supportHits());
	BOOST_CHECK_EQUAL(0, streaming.getBoundaryHandler().n_nonSupportHits());

	// TODO: finalize
	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_DataStructures_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryHandler_Couette_test) {
	pout << "SemiLagrangianBoundaryHandler_Couette_test..." << endl;

	double viscosity = 1e-1;
	double u0 = 0.1;
	size_t ref_level = 3;
	boost::shared_ptr<ProblemDescription<2> > couette = boost::make_shared<
			CouetteFlowGrad2D>(viscosity, u0, ref_level);

	// TODO: finalize
	pout << "done." << endl;
} /* SemiLagrangianBoundaryHandler_Couette_test */

BOOST_AUTO_TEST_SUITE_END()
