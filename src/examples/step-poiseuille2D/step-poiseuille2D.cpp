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
#include "natrium/utilities/CommandLineParser.h"
#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/benchmarks/PoiseuilleFlow2D.h"


using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-poiseuille2D..." << endl;
	CommandLineParser parser(argc, argv);



	parser.addDocumentationString("step-poiseuille2D",
			"2D Poiseuille Flow");
	parser.setPositionalArgument<int>("ref-level",
				"refinement of the computational grid");
	parser.setPositionalArgument<double>("ma","Mach number of the flow");
	parser.setArgument<double>("Re","Reynolds number of the flow",100.0);

	try {
		parser.importOptions();
	} catch (HelpMessageStop&) {
		return 0;
	}


	const double CFL = 0.4;
	const double Re = parser.getArgument<double>("Re");
	const double u_bulk = 0.0001 / 1.5; //1.0;
	const double height = 3.0;
	const double length = 2.0;
	const double orderOfFiniteElement = 2;
	//const double Ma = atof(argv[1]);
	//const double refinement_level = atoi(argv[2]);
	bool is_periodic = true;

	/// create CFD problem
	double viscosity  = u_bulk * height / Re;
	const double scaling = sqrt(3) * 1.5 * u_bulk / parser.getArgument<double>("ma");
	boost::shared_ptr<ProblemDescription<2> > poiseuille2D = boost::make_shared<
			PoiseuilleFlow2D>(viscosity,  parser.getArgument<int>("ref-level"), u_bulk, height,
			length, is_periodic);
	//viscosity = 0.5*dt*scaling*scaling/3.; //u_bulk * height / Re;
	//poiseuille2D->setViscosity(viscosity);
	//poiseuille2D->getExternalForce()->scale(viscosity);


	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/poiseuille2D";
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	//configuration->setRestartAtLastCheckpoint(false);

	configuration->setUserInteraction(false);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	configuration->setOutputTableInterval(10);
	configuration->setOutputCheckpointInterval(100000000);
	configuration->setOutputSolutionInterval(100);
	configuration->setCommandLineVerbosity(WELCOME);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setCFL(CFL);
	configuration->setForcingScheme(SHIFTING_VELOCITY);

	parser.applyToSolverConfiguration(*configuration);
	//configuration->setTimeIntegrator(OTHER);
	//configuration->setDealIntegrator(CRANK_NICOLSON);

	//configuration->setInitializationScheme(ITERATIVE);
	//configuration->setIterativeInitializationNumberOfIterations(100);
	//configuration->setIterativeInitializationResidual(1e-15);

	configuration->setConvergenceThreshold(1e-10);

	// make solver object and run simulation
	CFDSolver<2> solver(configuration, poiseuille2D);
	solver.run();

	pout << "Max Velocity  " <<
			solver.getMaxVelocityNorm() << "   (expected: "<< 1.5*u_bulk << ")" <<  endl;

	pout << "NATriuM step-poiseuille terminated." << endl;

	return 0;
}
