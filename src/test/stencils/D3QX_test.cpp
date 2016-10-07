/**
 * @file D3QX_test.cpp
 * @short
 * @date 11.03.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"

using std::exception;

using namespace natrium;

BOOST_AUTO_TEST_SUITE(D3QX_test)

BOOST_AUTO_TEST_CASE(D3QXMomentTrafo_test) {

	/////////////////
	// SANITY TEST //
	/////////////////

	vector<boost::shared_ptr<Stencil> > stencils;
	stencils.push_back(boost::make_shared<D3Q15>(1.0));
	stencils.push_back(boost::make_shared<D3Q15>(30.0));
	stencils.push_back(boost::make_shared<D3Q19>(1.0));
	stencils.push_back(boost::make_shared<D3Q19>(1000.0));
	stencils.push_back(boost::make_shared<D3Q27>(1.0));
	stencils.push_back(boost::make_shared<D3Q27>(1000.0));

	for (size_t s = 0; s < stencils.size(); s++) {
		// standard scaling
		const Stencil& dqmodel = *stencils.at(s);

		size_t Q = dqmodel.getQ();
		pout << "D3Q" << Q << "_Scaling" << dqmodel.getScaling()
				<< "_MomentTrafo_test..." << endl;

		numeric_matrix f_to_m(Q);
		numeric_matrix m_to_f(Q);
		dqmodel.getMomentBasis(f_to_m);
		dqmodel.getInverseMomentBasis(m_to_f);
		numeric_matrix Ident(Q);
		m_to_f.mmult(Ident, f_to_m);
		for (size_t i = 0; i < Q; i++) {
			for (size_t j = 0; j < Q; j++) {
				if (i != j)
					BOOST_CHECK_SMALL(Ident(i, j), 1e-10);
				else
					BOOST_CHECK_SMALL(Ident(i, j) - 1.0, 1e-10);
			}
		}

		pout << "done" << endl;
	}

} //D2Q9MomentTrafo_test

BOOST_AUTO_TEST_SUITE_END()

