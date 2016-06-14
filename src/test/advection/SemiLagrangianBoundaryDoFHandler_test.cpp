/*
 * SemiLagrangianBoundaryDoFHandler_test.cpp
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#include "natrium/advection/SemiLagrangianBoundaryDoFHandler.h"
#include "natrium/advection/SemiLagrangian.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"

#include "boost/test/unit_test.hpp"

namespace natrium {

BOOST_AUTO_TEST_SUITE(SemiLagrangianBoundaryDoFHandler_test)


BOOST_AUTO_TEST_CASE(BoundaryHit_Construction_test){
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

	BOOST_CHECK_NO_THROW(BoundaryHit<2>(dealii::Point<2>(0.0, 0.0), 0.0, dealii::Tensor<1, 2>(),
			LINEAR_RHO_U, cell, 1));

	pout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryDoFHandler_Construction_test){
	pout << "SemiLagrangianBoundaryDoFHandler_Construction_test..." << endl;

		BOOST_CHECK_NO_THROW(
				SemiLagrangianBoundaryDoFHandler<2> boundary_handler);

		pout << "done." << endl;
} /* SemiLagrangianBoundaryDoFHandler_Construction_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
