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
	cout << "PhysicalProperties_MassFluxX_test..." << endl;

	shared_ptr<ProblemDescription<2> > problem = make_shared<CouetteFlow2D> (1,1,2,1,40,true);
	shared_ptr<SolverConfiguration> config  = make_shared<SolverConfiguration>();
	config->setConvergenceThreshold(1e-4);
	config->setRestartAtLastCheckpoint(false);
	config->setUserInteraction(false);
	config->setSwitchOutputOff(true);
	config->setStencilScaling(1);
	config->setSedgOrderOfFiniteElement(2);
	const double dt = CFDSolverUtilities::calculateTimestep<2>(
			*problem->getTriangulation(), 2,
			D2Q9(1));
	config->setTimeStepSize(dt);

	CFDSolver<2> solver(config, problem);
	solver.run();

	const distributed_vector & ux = solver.getVelocity().at(0);
	double qx = PhysicalProperties<2>::meanVelocityX(ux, solver.getAdvectionOperator());
	BOOST_CHECK_CLOSE(qx, 0.5, 1e-4);

	cout << "done" << endl;
}


BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
