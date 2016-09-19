/**
 * @file BoundaryTools_test.cpp
 * @short Unit tests for all functions ins BoundaryTools.h
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/boundaries/BoundaryTools.h"

#include <string>

#include "boost/test/unit_test.hpp"

#include "deal.II/base/point.h"

#include "natrium/utilities/BasicNames.h"

namespace natrium {
namespace BoundaryTools{

BOOST_AUTO_TEST_SUITE(BoundaryTools_test)

BOOST_AUTO_TEST_CASE(BoundaryTools_CheckParallelLines_test) {

	pout << "BoundaryTools_CheckParallelLines_test..." << endl;

	dealii::Point<2> beginLine1(0.0, 0.0);
	dealii::Point<2> endLine1(0.0, 2.0);
	dealii::Point<2> beginLine2(1.0, 0.0);
	dealii::Point<2> endLine2(1.0, 2.0);
	std::string msg;

	// check if parallel lines could be detected
	BOOST_CHECK(checkParallelLines(beginLine1, endLine1, beginLine2, endLine2, msg));
	BOOST_CHECK(msg.empty());

	// check if non-parallel lines could be detected
	BOOST_CHECK(not checkParallelLines(beginLine1, endLine2, beginLine2, endLine1, msg));
    BOOST_CHECK(not msg.empty());

	// check if antiparallel lines are detected and respective line orientation is changed
	BOOST_CHECK(checkParallelLines(beginLine1, endLine1, endLine2, beginLine2, msg));
	BOOST_CHECK(beginLine2[1] == 2.0);
	BOOST_CHECK(endLine2[1] == 0.0);

	pout << "done." << endl;

} /* BoundaryTools_CheckParallelLines_test */

BOOST_AUTO_TEST_SUITE_END()
} /* namepace BoundaryTools */
} /* namespace natrium */
