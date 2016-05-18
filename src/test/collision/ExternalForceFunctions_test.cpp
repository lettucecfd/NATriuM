/*
 * ExternalForceFunctions_test.cpp
 *
 *  Created on: 12.05.2016
 *      Author: akraem3m
 */

#include <math.h>
#include <exception>

#include "boost/test/unit_test.hpp"

#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/collision/BGKStandard.h"
#include "natrium/utilities/Math.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/benchmarks/PoiseuilleFlow2D.h"
#include "natrium/benchmarks/PoiseuilleFlow3D.h"
#include "natrium/solver/CFDSolver.h"

using std::exception;

namespace natrium {

// is now an integration test

//BOOST_AUTO_TEST_SUITE(ExternalForceFunctions_test)
//
//BOOST_AUTO_TEST_CASE(ExternalForceFunctions_Poiseuille2D_test) {
//// this is rather an integration test than a Unit test
//	pout << "ExternalForceFunctions_Poiseuille2D_test..." << endl;
//
//	// setup test case
//	const double CFL = 0.8;
//	const double Re = 10;
//	const double u_bulk = 0.0001 / 1.5; //1.0;
//	const double height = 3.0;
//	const double length = 2.0;
//	const double orderOfFiniteElement = 2;
//	const double Ma = 0.1;
//	const double refinement_level = 1;
//	bool is_periodic = true;
//
//	/// create CFD problem
//	double viscosity = u_bulk * height / Re;
//	const double scaling = sqrt(3) * 1.5 * u_bulk / Ma;
//
//	/// setup configuration
//	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
//			SolverConfiguration>();
//	configuration->setSwitchOutputOff(true);
//	configuration->setUserInteraction(false);
//	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
//	configuration->setStencilScaling(scaling);
//	configuration->setCFL(CFL);
//	configuration->setConvergenceThreshold(1e-8);
//	//configuration->setTimeIntegrator(OTHER);
//	//configuration->setDealIntegrator(SDIRK_TWO_STAGES);
//
//	// SHIFTING VELOCITY
//	{
//		pout << "... shifting velocity scheme." << endl;
//		configuration->setForcingScheme(SHIFTING_VELOCITY);
//
//		// make solver object and run simulation
//		boost::shared_ptr<ProblemDescription<2> > poiseuille2D =
//				boost::make_shared<PoiseuilleFlow2D>(viscosity,
//						refinement_level, u_bulk, height, length, is_periodic);
//
//		CFDSolver<2> solver(configuration, poiseuille2D);
//		solver.run();
//
//		BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
//		pout << "... done." << endl;
//	}
//
//	// EXACT_DIFFERENCE
//	{
//		pout << "... exact difference scheme." << endl;
//		configuration->setForcingScheme(EXACT_DIFFERENCE);
//
//		// make solver object and run simulation
//		boost::shared_ptr<ProblemDescription<2> > poiseuille2D =
//				boost::make_shared<PoiseuilleFlow2D>(viscosity,
//						refinement_level, u_bulk, height, length, is_periodic);
//
//		CFDSolver<2> solver(configuration, poiseuille2D);
//		solver.run();
//
//		BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
//		pout << "... done." << endl;
//	}
//
//	//GUO
//	{
//		pout << "... Guo scheme." << endl;
//		configuration->setForcingScheme(GUO);
//
//		// make solver object and run simulation
//		boost::shared_ptr<ProblemDescription<2> > poiseuille2D =
//				boost::make_shared<PoiseuilleFlow2D>(viscosity,
//						refinement_level, u_bulk, height, length, is_periodic);
//		CFDSolver<2> solver(configuration, poiseuille2D);
//		solver.run();
//
//		BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
//		pout << "... done." << endl;
//	}
//
//	pout << "done" << endl;
//} /* ExternalForceFunctions_Poiseuille2D_test*/
//
//BOOST_AUTO_TEST_CASE(ExternalForceFunctions_Poiseuille3D_test) {
//// this is rather an integration test than a Unit test
//	pout << "ExternalForceFunctions_Poiseuille3D_test..." << endl;
//
//	// setup test case
//	const double CFL = 0.8;
//	const double Re = 10;
//	const double u_bulk = 0.0001 / 1.5; //1.0;
//	const double height = 3.0;
//	const double width = 2.0;
//	const double length = 2.0;
//	const double orderOfFiniteElement = 2;
//	const double Ma = 0.1;
//	const double refinement_level = 1;
//	bool is_periodic = true;
//
//	/// create CFD problem
//	double viscosity = u_bulk * height / Re;
//	const double scaling = sqrt(3) * 1.5 * u_bulk / Ma;
//
//	/// setup configuration
//	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
//			SolverConfiguration>();
//	configuration->setSwitchOutputOff(true);
//	configuration->setUserInteraction(false);
//	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
//	configuration->setStencilScaling(scaling);
//	configuration->setCFL(CFL);
//	configuration->setConvergenceThreshold(1e-8);
//	//configuration->setTimeIntegrator(OTHER);
//	//configuration->setDealIntegrator(SDIRK_TWO_STAGES);
//
//	// forall stencil types and forcing schemes
//	// TODO there is a bug in D3Q27; set "i < 4" when resolved
//	for (size_t i = 1; i < 3; i++) {
//		configuration->setStencil(static_cast<StencilType>(i));
//		for (size_t j = 1; j < 4; j++) {
//			configuration->setForcingScheme(static_cast<ForceType>(j));
//			switch (i) {
//			case 1: {
//				pout << "...D3Q19; ";
//				break;
//			}
//			case 2: {
//				pout << "...D3Q15; ";
//				break;
//			}
//			case 3: {
//				pout << "...D3Q27; ";
//				break;
//			}
//			}
//			switch (j) {
//			case 1: {
//				pout << "...shifting velocity scheme; ";
//				break;
//			}
//			case 2: {
//				pout << "...exact difference scheme; ";
//				break;
//			}
//			case 3: {
//				pout << "...Guo scheme; ";
//				break;
//			}
//			}
//
//			// make solver object and run simulation
//			boost::shared_ptr<ProblemDescription<3> > poiseuille3D =
//					boost::make_shared<PoiseuilleFlow3D>(viscosity,
//							refinement_level, u_bulk, height, width, length,
//							is_periodic);
//
//			CFDSolver<3> solver(configuration, poiseuille3D);
//			solver.run();
//
//			BOOST_CHECK_CLOSE(solver.getMaxVelocityNorm(), 1.5 * u_bulk, 5);
//			pout << "...... done." << endl;
//		}
//	}
//	pout << "done" << endl;
//} /* ExternalForceFunctions_Poiseuille3D_test*/
//
//BOOST_AUTO_TEST_SUITE_END()

}
