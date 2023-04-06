/**
 * @file D2Q25H_test.cpp
 * @short
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/stencils/D2Q19H.h"

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"

using std::exception;

using namespace natrium;

BOOST_AUTO_TEST_SUITE(D2Q19H_test)

    BOOST_AUTO_TEST_CASE(D2Q19HConstruction_test) {
        pout << "D2Q19HConstruction_test..." << endl;

        /////////////////
        // SANITY TEST //
        /////////////////
        D2Q19H dqmodel;
        BOOST_CHECK_EQUAL(dqmodel.getD(), size_t(2));
        BOOST_CHECK_EQUAL(dqmodel.getQ(), size_t(19));
        BOOST_CHECK_CLOSE(dqmodel.getSpeedOfSound(), sqrt(1./3.),0.000001);
        BOOST_CHECK_CLOSE(dqmodel.getSpeedOfSoundSquare(), 1./3.,0.000001);
        std::vector<double> weights = dqmodel.getWeights();
        double sum =0.0;
        for(int i = 0; i < dqmodel.getQ(); i++)
            sum+=weights[i];
        BOOST_CHECK_CLOSE(sum,1.0,0.0000001);

        pout << "done" << endl;
    } //D2Q19HConstruction_test



//////////////////////////////////
// TESTS FOR THE SCALED VERSION //
//////////////////////////////////
    const double SCALING = 0.5;

    BOOST_AUTO_TEST_CASE(D2Q19HConstruction_Scaled_test) {
        pout << "D2Q19HConstruction_Scaled_test..." << endl;

        /////////////////
        // SANITY TEST //
        /////////////////
        D2Q19H dqmodel(SCALING);
        BOOST_CHECK_EQUAL(dqmodel.getD(), size_t(2));
        BOOST_CHECK_EQUAL(dqmodel.getQ(), size_t(19));
        BOOST_CHECK_CLOSE(dqmodel.getSpeedOfSound(), SCALING*sqrt(1./3.),0.000001);
        BOOST_CHECK_CLOSE(dqmodel.getSpeedOfSoundSquare(), SCALING*SCALING*1./3.,0.000001);

        pout << "done" << endl;
    } //D2Q19HConstruction_Scaled_test



BOOST_AUTO_TEST_SUITE_END()


