/**
 * @file DirichletBoundaryU2D_test.cpp
 * @short
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <natrium/problemdescription/DirichletBoundaryU.h>
#include "boost/test/unit_test.hpp"

#include "deal.II/base/function.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/fe/component_mask.h"
#include "deal.II/base/index_set.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/advection/SEDGMinLee.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "WallTestDomain2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(DirichletBoundaryU2D_test)


BOOST_AUTO_TEST_CASE(DirichletBoundaryU2D_MassConservation_test) {
	pout << "DirichletBoundaryRhoU2D_MassConservation_test..." << endl;

	// make problem and solver
	boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<WallTestDomain2D>(
			1, true);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("/tmp");
	configuration->setSwitchOutputOff(true);
	configuration->setNumberOfTimeSteps(100);
	configuration->setSedgOrderOfFiniteElement(1);
	configuration->setTimeStepSize(0.01);

	CFDSolver<2> solver(configuration, problem);

	solver.run();

	// check mass conservation
	double mass = solver.getDensity().l1_norm();
	mass /= solver.getNumberOfDoFs();
	BOOST_CHECK_SMALL(mass - 1.0, 1e-10);

	pout << "done" << endl;
} /*DirichletBoundaryU2D_MassConservation_test */

BOOST_AUTO_TEST_CASE(DirichletBoundaryU2D_CompareToRhoU_test) {
	pout << "DirichletBoundaryU2D_CompareToRhoU_test..." << endl;
	// If we perform only one explicit Euler step the two solutions have to agree
	// construct test problem with prescribed rho and u
	boost::shared_ptr<ProblemDescription<2> > problem1 = boost::make_shared<WallTestDomain2D>(
			1);
	// construct test problem with prescribed u
	boost::shared_ptr<ProblemDescription<2> > problem2 = boost::make_shared<WallTestDomain2D>(
				1, true);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(true);
	configuration->setUserInteraction(false);
	configuration->setNumberOfTimeSteps(1);
	configuration->setSedgOrderOfFiniteElement(1);
	configuration->setTimeIntegrator(OTHER);
	configuration->setDealIntegrator(FORWARD_EULER);
	configuration->setTimeStepSize(0.1);

	CFDSolver<2> solver1(configuration, problem1);

	std::ofstream fout("/tmp/natrium_matrices.txt");
	solver1.getAdvectionOperator()->getSystemMatrix().print(fout);
	solver1.getAdvectionOperator()->getSystemVector().print(fout);
	solver1.stream();

	CFDSolver<2> solver2(configuration, problem2);

	solver2.getAdvectionOperator()->getSystemMatrix().print(fout);
	solver2.stream();
	fout.close();

	// check equality
	for (size_t i = 0; i < solver1.getNumberOfDoFs(); i++){
		for (size_t j = 0; j < solver1.getStencil()->getQ(); j++){
			cout << "i,j:" << i << " " << j << endl;
			BOOST_CHECK_SMALL(solver1.getF().at(j)(i) - solver2.getF().at(j)(i), 1e-10);
		}
	}

	pout << "done" << endl;
} /*DirichletBoundaryU2D_CompareToRhoU_test */

BOOST_AUTO_TEST_SUITE_END() /*DirichletBoundaryU2D_test*/

} /* namespace natrium */
