/*
 * ExternalForceFunctions_test.cpp
 *
 *  Created on: 12.05.2016
 *      Author: akraem3m
 */



#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/collision/BGKStandard.h"
#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/benchmarks/PoiseuilleFlow2D.h"

using std::exception;

namespace natrium {

BOOST_AUTO_TEST_SUITE(ExternalForceFunctions_test)

BOOST_AUTO_TEST_CASE(ExternalForceFunctions_Poiseuille2D_test) {

	pout << "ExternalForceFunctions_Poiseuille2D_test..." << endl;


	pout << "done" << endl;
} /* ExternalForceFunctions_Poiseuille2D_test*/

BOOST_AUTO_TEST_SUITE_END()

}
