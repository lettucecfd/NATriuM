#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q13.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q21.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/stencils/RD3Q27.h"
#include "natrium/stencils/D2Q25H.h"
#include "natrium/stencils/D2Q19V.h"
#include "natrium/stencils/D2Q19H.h"
#include "natrium/stencils/D3Q45.h"
#include "natrium/stencils/D3Q77.h"
#include "natrium/stencils/D3V27.h"
#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"

using std::exception;

using namespace natrium;

template<class T_stencil>
void stencil_test_function()
{
    T_stencil stencil;
    pout << "D" << stencil.getD() << "Q" << stencil.getQ() << "_Stencil_test...\n" ;
    BOOST_CHECK_CLOSE(stencil.getSpeedOfSound(), sqrt(1./3.),0.000001);
    BOOST_CHECK_CLOSE(stencil.getSpeedOfSoundSquare(), 1./3.,0.000001);

    // Check weight sum equals zero
    const std::vector<double> weights = stencil.getWeights();
    double sum = 0.0;
    for(size_t i = 0; i < stencil.getQ(); i++)
        sum+=weights[i];
    BOOST_CHECK_SMALL(1.0-sum,0.0000001);

    // Check velocity sum equals zero in every direction
    for (size_t dims = 0; dims < stencil.getD(); dims++) {
        double sum_speeds = 0.0;
        for(size_t i = 0; i < stencil.getQ(); i++) {
            sum_speeds += stencil.getDirection(i)[dims];
        }
        BOOST_CHECK_SMALL(0.0-sum_speeds, 0.0000001);
    }
    pout << "done.\n" ;
}


BOOST_AUTO_TEST_SUITE(Stencil_test_suite)

BOOST_AUTO_TEST_CASE(Stencil_tests) {
        stencil_test_function<D2Q9>();
        stencil_test_function<D2Q19H>();
        stencil_test_function<D2Q19V>();
        stencil_test_function<D2Q25H>();
        stencil_test_function<D3Q13>();
        stencil_test_function<D3Q15>();
        stencil_test_function<D3Q19>();
        stencil_test_function<D3Q21>();
        stencil_test_function<D3Q27>();
        stencil_test_function<D3V27>();
        stencil_test_function<D3Q45>();
        stencil_test_function<D3Q77>();
    }

BOOST_AUTO_TEST_SUITE_END()
