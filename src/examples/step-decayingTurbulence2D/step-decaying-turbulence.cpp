/**
 * @file step-2.cpp
 * @short Second tutorial:  Poiseuille Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"

#include "DecayingTurbulence2D.h"
#include "natrium/utilities/CFDSolverUtilities.h"


using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-decaying-turbulence..." << endl;

	const double CFL = 0.4;
	const double Re = 10000;
	const double u_bulk = 0.0001 / 1.5; //1.0;
	const double orderOfFiniteElement = 2;
	const double Ma = atof(argv[1]);
	const double refinement_level = atoi(argv[2]);

	/// create CFD problem
	double viscosity  = u_bulk * 1.0 / Re;
	const double scaling = sqrt(3) * 1.5 * u_bulk / Ma;
	boost::shared_ptr<ProblemDescription<2> > decaying = boost::make_shared<
			DecayingTurbulence2D>(viscosity, refinement_level);
	//viscosity = 0.5*dt*scaling*scaling/3.; //u_bulk * height / Re;
	//poiseuille2D->setViscosity(viscosity);
	//poiseuille2D->getExternalForce()->scale(viscosity);


	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/decaying-turbulence2D";
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(10);
	configuration->setOutputCheckpointInterval(100000000);
	configuration->setOutputSolutionInterval(100);
	configuration->setCommandLineVerbosity(WELCOME);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setCFL(CFL);
	//configuration->setTimeIntegrator(OTHER);
	//configuration->setDealIntegrator(CRANK_NICOLSON);

	//configuration->setInitializationScheme(ITERATIVE);
	//configuration->setIterativeInitializationNumberOfIterations(100);
	//configuration->setIterativeInitializationResidual(1e-15);

	configuration->setConvergenceThreshold(1e-10);

	// make solver object and run simulation
	CFDSolver<2> solver(configuration, decaying);
	solver.run();

	pout << "Max Velocity  " <<
			solver.getMaxVelocityNorm() << "   (laminar: "<< 1.5*u_bulk << ")" <<  endl;

	pout << "NATriuM step-turbulent-channel terminated." << endl;

	return 0;
}
