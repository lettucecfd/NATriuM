/**
 * @file MinLeeBoundary2D_test.cpp
 * @short
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <problemdescription/MinLeeBoundary.h>

#include <boost/test/unit_test.hpp>

#include "deal.II/base/function.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/fe/component_mask.h"

#include "utilities/BasicNames.h"
#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "problemdescription/ProblemDescription.h"
#include "advection/SEDGMinLee.h"
#include "WallTestDomain2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(MinLeeBoundary2D_test)

class BoundaryTestDensity: public dealii::Function<2> {
public:
	virtual double value(const dealii::Point<2> &p,
			const unsigned int component = 0) const {
		return 1;
	}
};
class BoundaryTestVelocity: public dealii::Function<2> {
public:
	virtual void vector_value(const dealii::Point<2> &p,
			dealii::Vector<double> &values) const {
		values(0) = 0;
		values(1) = 0;
	}
};

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_Construction_test) {
	cout << "MinLeeBoundary2D_Construction_test..." << endl;

	BOOST_CHECK_NO_THROW(MinLeeBoundary<2> mlBound1(0, make_shared<BoundaryTestDensity>(), make_shared<BoundaryTestVelocity>()));
	numeric_vector U(2);
	BOOST_CHECK_NO_THROW(MinLeeBoundary<2> mlBound2(0, U);
);

	cout << "done" << endl;
} /*MinLeeBoundary2D_Construction_test */

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_SparsityPattern_test){
	cout << "MinLeeBoundary2D_SparsityPattern_test..." << endl;

	shared_ptr<ProblemDescription<2> > problem = make_shared<WallTestDomain2D>(1);
	SEDGMinLee<2> advector(problem->getTriangulation(),
			problem->getBoundaries(),
			2, make_shared<D2Q9IncompressibleModel>());
	vector< bool > isBoundary;
	for (size_t i = 0; i++; i < advector.getNumberOfDoFs()){
		// left boundary
		std::set<dealii::types::boundary_id> boundaryIndicators;
		boundaryIndicators.insert(0);
		dealii::DoFTools::extract_boundary_dofs(*(advector.getDoFHandler()), dealii::ComponentMask(), isBoundary, boundaryIndicators);
		if (isBoundary.at(i)){
			BOOST_CHECK(advector.getBlockSparsityPattern().block(4,6).exists(i,i));
			BOOST_CHECK(advector.getBlockSparsityPattern().block(0,2).exists(i,i));
			BOOST_CHECK(advector.getBlockSparsityPattern().block(5,7).exists(i,i));
		}
	}

	cout << "done" << endl;
} /* MinLeeBoundary2D_SparsityPattern_test */

BOOST_AUTO_TEST_CASE(MinLeeBoundary2D_MassConversion_test) {
	cout << "MinLeeBoundary2D_MassConversion_test..." << endl;

	cout << "done" << endl;
} /*MinLeeBoundary2D_MassConversion_test */

BOOST_AUTO_TEST_SUITE_END() /*MinLeeBoundary2D_test*/

} /* namespace natrium */
