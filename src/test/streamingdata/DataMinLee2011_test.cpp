/**
 * @file DataMinLee2011_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "streamingdata/DataMinLee2011.h"

#include "boost/test/unit_test.hpp"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "PeriodicFlow2D.h"

#include "utilities/BasicNames.h"


namespace natrium {

BOOST_AUTO_TEST_SUITE(DataMinLee2011_test)

BOOST_AUTO_TEST_CASE(DataMinLee2011_Construction_test) {
	cout << "DataMinLee2011_Construction_test..." << endl;

	double relaxationParameter = 0.7;
	numeric_vector velocity(2);
	velocity(0) = 0.05;
	velocity(1) = 0.01;
	size_t fe_order = 2;
	PeriodicFlow2D periodic(relaxationParameter, velocity);
	DataMinLee2011<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>());

	cout << "done.";
} /* DataMinLee2011_Construction_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
