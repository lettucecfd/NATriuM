/**
 * @file DataMinLee2011_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "streamingdata/DataMinLee2011.h"

#include "boost/test/unit_test.hpp"

#include "PeriodicFlow2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(DataMinLee2011_test)

BOOST_AUTO_TEST_CASE(DataMinLee2011_Construction_test) {
	cout << "DataMinLee2011_Construction_test..." << endl;


	/*DataMinLee2011(shared_ptr<dealii::Triangulation<dim> > triangulation,
			shared_ptr<BoundaryCollection<dim> > boundaries,
			size_t orderOfFiniteElement,
			shared_ptr<BoltzmannModel> boltzmannModel);
*/
	cout << "done.";
} /* DataMinLee2011_Construction_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
