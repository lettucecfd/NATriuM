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
        BOOST_CHECK_EQUAL(dqmodel.D, size_t(2));
        BOOST_CHECK_EQUAL(dqmodel.Q, size_t(19));
        BOOST_CHECK_CLOSE(dqmodel.getSpeedOfSound(), sqrt(1./3.),0.000001);
        BOOST_CHECK_CLOSE(dqmodel.getSpeedOfSoundSquare(), 1./3.,0.000001);

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



    BOOST_AUTO_TEST_CASE(D2Q19HMomentTrafo_test) {
            pout << "D2Q19HMomentTrafo_test..." << endl;

            /////////////////
            // SANITY TEST //
            /////////////////

            // standard scaling
            D2Q19H dqmodel(1.0);

            numeric_matrix f_to_m(19);
            numeric_matrix m_to_f(19);
            dqmodel.getMomentBasis(f_to_m);
            dqmodel.getInverseMomentBasis(m_to_f);
            numeric_matrix Ident(19);
            m_to_f.mmult(Ident, f_to_m);
            for (size_t i = 0; i < 19; i++){
                for (size_t j = 0; j < 19; j++) {
                    if (i != j)
                        BOOST_CHECK_SMALL(Ident(i, j), 1e-10);
                    else
                        BOOST_CHECK_SMALL(Ident(i, j) - 1.0, 1e-10);
                }
            }

            D2Q19H dqmodel2(1000.0);

            dqmodel.getMomentBasis(f_to_m);
            dqmodel.getInverseMomentBasis(m_to_f);
            m_to_f.mmult(Ident, f_to_m);
            for (size_t i = 0; i < 19; i++){
                for (size_t j = 0; j < 19; j++) {
                    if (i != j)
                        BOOST_CHECK_SMALL(Ident(i, j), 1e-10);
                    else
                        BOOST_CHECK_SMALL(Ident(i, j) - 1.0, 1e-10);
                }
            }


            pout << "done" << endl;

    }
BOOST_AUTO_TEST_SUITE_END()


