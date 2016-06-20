/*
 * SemiLagrangianBoundaryDoFHandler_test.cpp
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#include "natrium/advection/SemiLagrangianBoundaryDoFHandler.h"

#include "boost/test/unit_test.hpp"

namespace natrium {

BOOST_AUTO_TEST_SUITE(SemiLagrangianBoundaryDoFHandler_test)

BOOST_AUTO_TEST_CASE(SemiLagrangianBoundaryDoFHandler_Construction_test){
	pout << "SemiLagrangianBoundaryDoFHandler_Construction_test..." << endl;

		BOOST_CHECK_NO_THROW(
				SemiLagrangianBoundaryDoFHandler<2> boundary_handler);

		pout << "done." << endl;
} /* SemiLagrangianBoundaryDoFHandler_Construction_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
