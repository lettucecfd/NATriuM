/*
 * SemiLagrangianBoundaryDoFHandler_test.cpp
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#include "natrium/advection/SemiLagrangianBoundaryDoFHandler.h"
#include "natrium/advection/SemiLagrangian.h"
#include "natrium/problemdescription/BoundaryCollection.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/problemdescription/LinearBoundaryRhoU.h"
#include "natrium/advection/SemiLagrangianVectorReferenceTypes.h"

#include "boost/test/unit_test.hpp"

namespace natrium {

BOOST_AUTO_TEST_SUITE(SemiLagrangianBoundaryDoFHandler_test)

BOOST_AUTO_TEST_CASE(BoundaryHit_Construction_test) {
	pout << "BoundaryHit_Construction_test" << endl;

	// make dof handler
	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();

	OutgoingDistributionValue a(false, 0, 1);
	//BOOST_CHECK_NO_THROW(
	BoundaryHit<2>(dealii::Point<2>(0.0, 0.0), 0.0, dealii::Tensor<1, 2>(),
			*(periodic.getBoundaries()->getBoundary(0)), cell, a);
	//);

	pout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryDoFHandler_Construction_test) {
	pout << "SemiLagrangianBoundaryDoFHandler_Construction_test..." << endl;

	// make dof handler

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	dealii::DoFHandler<2> dof_handler(*periodic.getMesh());
	dealii::FE_DGQArbitraryNodes<2> fe(dealii::QGaussLobatto<1>(fe_order + 1));
	dof_handler.distribute_dofs(fe);

	D2Q9 d2q9;
	SemiLagrangianBoundaryDoFHandler<2>(dof_handler.locally_owned_dofs(), d2q9);

	dof_handler.clear();
	pout << "done." << endl;
} /* SemiLagrangianBoundaryDoFHandler_Construction_test */

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryDoFHandler_Apply_test) {
	pout << "SemiLagrangianBoundaryDoFHandler_Apply_test..." << endl;

	// make dof handler

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	numeric_vector u(2);
	u[0] = 0.1;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
	LinearBoundaryRhoU<2> wall(0, u);
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();

	D2Q9 d2q9;
	SemiLagrangianBoundaryDoFHandler<2> sl_handler(
			sl.getDoFHandler()->locally_owned_dofs(), d2q9);

	vector<distributed_vector> v1;
	for (size_t i = 0; i < 9; i++) {
		distributed_vector tmp(sl.getDoFHandler()->locally_owned_dofs(),
				MPI_COMM_WORLD);
		v1.push_back(tmp);
	}
	DistributionFunctions fold(v1);

	vector<distributed_vector> v2;
	for (size_t i = 0; i < 9; i++) {
		distributed_vector tmp(sl.getDoFHandler()->locally_owned_dofs(),
				MPI_COMM_WORLD);
		v2.push_back(tmp);
	}
	DistributionFunctions fnew(v2);

	//////////////////////////////////////////////////////////
	//    b.hit2    2ndary_dof       b.hit1        dof		//
	//    routine        data        routine       data     //
	// ......----->        X       .....----->      X		//
	//  (5)   (2)         (3)      (4)    (1)				//
	//////////////////////////////////////////////////////////
	// add boundary hit (1)
	OutgoingDistributionValue a(false,
			sl.getDoFHandler()->locally_owned_dofs().nth_index_in_set(0), 1);
	dealii::Tensor<1, 2> n;
	n[0] = 1;
	BoundaryHit<2> hit1(dealii::Point<2>(0.0, 0.0), 0.0, n, wall, cell, a);
	LagrangianPathDestination destination1 = sl_handler.addBoundaryHit(hit1);

	// add a secondary boundary hit
	OutgoingDistributionValue a2;
	// (2)
	BoundaryHit<2> hit2(dealii::Point<2>(0.0, 0.0), 0.0, n, wall, cell, a2);
	// (3)
	LagrangianPathDestination destination2 = sl_handler.addBoundaryHit(hit2);
	// configure secondary as input to primary (4)
	BOOST_CHECK(destination2.isBoundaryHit);
	IncomingDistributionValue fv(
			sl_handler.getBoundaryHit(destination2).out.getIndex());
	sl_handler.getBoundaryHit(destination1).in.push_back(fv);

	// configure input of secondary boundary hit (5)
	IncomingDistributionValue fv2;
	fv2.internalDoFs.push_back(
			sl.getDoFHandler()->locally_owned_dofs().nth_index_in_set(0));
	fv2.internalValues.push_back(1.0);
	sl_handler.getBoundaryHit(destination2).in.push_back(fv2);

	fold.at(1)(sl.getDoFHandler()->locally_owned_dofs().nth_index_in_set(0)) = 1;
	sl_handler.calculateAndApplyBoundaryValues(fnew, fold, 0.1);


	// calculate expected
	/*f[boundary_hit.out] = f(boundary_hit.in.at(0))
			+ 2.0 * stencil.getWeight(boundary_hit.out.getAlpha()) * 1.0
					* (ea * velocity) / stencil.getSpeedOfSoundSquare();
	  	  f_secondary_bv = 1.0 + 2.0 * 1/9 * 1.0 * u0 / ( 1/3) = 32/30
	  	  f_out          = 32/30 + 2.0 * 1/9 * 1.0 * u0 / ( 1/3) = 34/30
	*/
	//BOOST_CHECK_SMALL(fnew.at(1)(sl.getDoFHandler()->locally_owned_dofs().nth_index_in_set(0)) - 34.0/30.0, 1e-10);

	pout << "done." << endl;
} /* SemiLagrangianBoundaryDoFHandler_Apply_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
