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
#include "natrium/benchmarks/ShearLayer2D.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/solver/CFDSolver.h"

using namespace natrium;

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

BOOST_AUTO_TEST_CASE(PhysicalProperties_Mass_test) {
	pout << "PhysicalProperties_Mass_test..." << endl;

	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<CouetteFlow2D> (1,1,2,1,40,true);
	boost::shared_ptr<SolverConfiguration> config  = boost::make_shared<SolverConfiguration>();
	config->setConvergenceThreshold(1e-4);
	config->setUserInteraction(false);
	config->setSwitchOutputOff(true);
	config->setStencilScaling(1);
	config->setSedgOrderOfFiniteElement(2);

	config->setCFL(0.4);

	CFDSolver<2> solver(config, problem);

	double m = PhysicalProperties<2>::mass(solver.getDensity(), solver.getAdvectionOperator());
	BOOST_CHECK_CLOSE(m, 1, 1e-4);

	pout << "done" << endl;
} /* PhysicalProperties_Mass_test*/

BOOST_AUTO_TEST_CASE(PhysicalProperties_KineticEnergy_test) {
	pout << "PhysicalProperties_KineticEnergy_test..." << endl;

	// Setup as in Minion and Brown's paper from 1995
	const double u0 = 1;
	const double kappa = 80;
	const double viscosity = 0.0001;
	const size_t refinement_level = 3;

	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<ShearLayer2D> (viscosity, refinement_level, u0, kappa);
	boost::shared_ptr<SolverConfiguration> config  = boost::make_shared<SolverConfiguration>();
	config->setUserInteraction(false);
	config->setSwitchOutputOff(true);
	config->setStencilScaling(10);
	config->setSedgOrderOfFiniteElement(2);

	CFDSolver<2> solver(config, problem);

	const vector<distributed_vector> & u = solver.getVelocity();
	const distributed_vector & rho = solver.getDensity();
	double E = PhysicalProperties<2>::kineticEnergy(u,rho, solver.getAdvectionOperator());
	BOOST_CHECK_CLOSE(E, 0.476, 10);

	pout << "done" << endl;
} /* PhysicalProperties_KineticEnergy_test */

BOOST_AUTO_TEST_CASE(PhysicalProperties_Enstrophy_test) {
	pout << "PhysicalProperties_Enstrophy_test..." << endl;

	// Setup as in Minion and Brown's paper from 1995
	const double u0 = 1;
	const double kappa = 80;
	const double viscosity = 0.0001;
	const size_t refinement_level = 4;

	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<ShearLayer2D> (viscosity, refinement_level, u0, kappa);
	boost::shared_ptr<SolverConfiguration> config  = boost::make_shared<SolverConfiguration>();
	config->setUserInteraction(false);
	config->setSwitchOutputOff(true);
	config->setStencilScaling(10);
	config->setSedgOrderOfFiniteElement(9);

	CFDSolver<2> solver(config, problem);

	const vector<distributed_vector> & u = solver.getVelocity();
	double sq;
	double E = PhysicalProperties<2>::enstrophy(u, solver.getAdvectionOperator(), &sq);
	BOOST_CHECK_CLOSE(E, 213, 10);

	pout << "done" << endl;
} /* PhysicalProperties_Enstrophy_test */


BOOST_AUTO_TEST_CASE(PhysicalProperties_EntropyGhosted_test) {
	pout << "PhysicalProperties_EntropyGhosted_test..." << endl;

	// Setup as in Minion and Brown's paper from 1995
	const double u0 = 1;
	const double kappa = 80;
	const double viscosity = 0.0001;
	const size_t refinement_level = 7;

	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<ShearLayer2D> (viscosity, refinement_level, u0, kappa);
	boost::shared_ptr<SolverConfiguration> config  = boost::make_shared<SolverConfiguration>();
	config->setUserInteraction(false);
	config->setSwitchOutputOff(true);
	config->setStencilScaling(1);
	config->setSedgOrderOfFiniteElement(1);
	config->setAdvectionScheme(SEMI_LAGRANGIAN);

	CFDSolver<2> solver(config, problem);

	const DistributionFunctions& f = solver.getF();
	BOOST_CHECK_NO_THROW(PhysicalProperties<2>::entropy(f, solver.getAdvectionOperator()));



	pout << "done" << endl;
} /* PhysicalProperties_EntropyGhosted_test */



BOOST_AUTO_TEST_SUITE_END()

