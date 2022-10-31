/**
 * @file step-2.cpp
 * @short Second tutorial:  Poiseuille Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/D3Q45.h"
#include "natrium/utilities/CommandLineParser.h"
#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "BoundaryTemperature.h"

#include "natrium/benchmarks/PoiseuilleFlow3D.h"


using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-poiseuille3D..." << endl;
	CommandLineParser parser(argc, argv);



	parser.addDocumentationString("step-poiseuille3D",
			"3D Poiseuille Flow");
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
	const double u_bulk = 1.0;
	const double height = 1.0;
    const double width = 1.0;
	const double length = 1.0;
	const double orderOfFiniteElement = 2;
	//const double Ma = atof(argv[1]);
	//const double refinement_level = atoi(argv[2]);
	bool is_periodic = true;

	/// create CFD problem
	double viscosity  = u_bulk * height / Re;
	const double scaling = sqrt(3) * u_bulk / parser.getArgument<double>("ma");
	boost::shared_ptr<ProblemDescription<3> > poiseuille3D = boost::make_shared<
			PoiseuilleFlow3D>(viscosity,  parser.getArgument<int>("ref-level"), u_bulk, height, width,
			length, is_periodic);
	//viscosity = 0.5*dt*scaling*scaling/3.; //u_bulk * height / Re;
	//poiseuille3D->setViscosity(viscosity);
	//poiseuille3D->getExternalForce()->scale(viscosity);


	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/poiseuille3D";
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
    configuration->setHeatCapacityRatioGamma(1.4);


    parser.applyToSolverConfiguration(*configuration);
	//configuration->setTimeIntegrator(OTHER);
	//configuration->setDealIntegrator(CRANK_NICOLSON);

	//configuration->setInitializationScheme(ITERATIVE);
	//configuration->setIterativeInitializationNumberOfIterations(100);
	//configuration->setIterativeInitializationResidual(1e-15);

	configuration->setConvergenceThreshold(1e-10);

	// make solver object and run simulation

	CompressibleCFDSolver<3> solver(configuration, poiseuille3D);
    solver.appendDataProcessor(
            boost::make_shared<BoundaryTemperature>(solver));
    solver.run();

	pout << "Max Velocity  " <<
			solver.getMaxVelocityNorm() << "   (expected: "<< 1.5*u_bulk << ")" <<  endl;

	pout << "NATriuM step-poiseuille terminated." << endl;

	return 0;
}
