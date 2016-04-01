/*
 * PhysicalProperties_test.cpp
 *
 *  Created on: Jun 6, 2014
 *      Author: kraemer
 */

#include "natrium/solver/PhysicalProperties.h"


#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"
#include "natrium/benchmarks/CouetteFlow2D.h"
#include "natrium/utilities/CFDSolverUtilities.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(PhysicalProperties_test)

BOOST_AUTO_TEST_CASE(PhysicalProperties_MassFluxX_test) {
	pout << "PhysicalProperties_MassFluxX_test..." << endl;

	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<CouetteFlow2D> (1,1,2,1,40,true);
	boost::shared_ptr<SolverConfiguration> config  = boost::make_shared<SolverConfiguration>();
	config->setConvergenceThreshold(1e-4);
	config->setUserInteraction(false);
	config->setSwitchOutputOff(true);
	config->setStencilScaling(1);
	config->setSedgOrderOfFiniteElement(2);

	config->setCFL(0.4);

	CFDSolver<2> solver(config, problem);
	solver.run();

	const distributed_vector & ux = solver.getVelocity().at(0);
	double qx = PhysicalProperties<2>::meanVelocityX(ux, solver.getAdvectionOperator());
	BOOST_CHECK_CLOSE(qx, 0.5, 1e-4);

	pout << "done" << endl;
}


BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
