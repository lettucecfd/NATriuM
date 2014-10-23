/**
 * @file Collision_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "collision/Collision.h"


#include "collision/BGKTransformed.h"

#include "boost/test/unit_test.hpp"

#include "solver/DistributionFunctions.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "boltzmannmodels/BoltzmannModel.h"

#include "utilities/BasicNames.h"


namespace natrium {

BOOST_AUTO_TEST_SUITE(Collision_test)


BOOST_AUTO_TEST_CASE(Collision_test_low_level_D2Q9) {
/*
	BGKTransformed bgk(9, 0.5);
	D2Q9IncompressibleModel d2q9(0.5);

	DistributionFunctions f;
	f.reinit(9, 5);


	d2q9.getEquilibriumDistribution(i,u,rho);
	bgk.collideSingleDoF()
	Dis*/


}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
