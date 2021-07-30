/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>
#include <mpi.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/CompressibleCFDSolver.h"

#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/CommandLineParser.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/stencils/D2Q9.h"

#include "natrium/utilities/BasicNames.h"

#include "Cylinder2D.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

    MPIGuard::getInstance();
    CommandLineParser parser(argc, argv);
    parser.setArgument<double>("Ma", "Mach number", 0.1);
    parser.setArgument<int>("Re", "Reynolds number", 100);

    parser.setArgument<int>("compressible", "Compressible CFD solver needed?", 0);

    try {
        parser.importOptions();
    } catch (HelpMessageStop&){
        return 0;
    }


	pout << "Starting NATriuM step-9..." << endl;


	// set Reynolds and Mach number
	const double Ma = parser.getArgument<double>("Ma");
    const double Re = parser.getArgument<int>("Re");;


    // set Problem so that the right Re and Ma are achieved
	const double U = 1/sqrt(3)*Ma;
	const double dqScaling = 1;
	const double viscosity = U / Re; // (because L = 1)


	// load grid
	boost::shared_ptr<Cylinder2D> cylinder = boost::make_shared<Cylinder2D>(
				viscosity, U);
	D2Q9 stencil(dqScaling);
	// set FE order and time step size
	const size_t orderOfFiniteElement = 2;
	const double cfl=5;


	pout << "Mach number: " << U / ( dqScaling / sqrt(3)) << endl;

	// configure solver
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-9";
	configuration->setOutputDirectory(dirname.str());
	configuration->setUserInteraction(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setSimulationEndTime(5000);
	configuration->setTimeIntegrator(EXPONENTIAL);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setCFL(cfl);
	configuration->setCommandLineVerbosity(7);
	configuration->setHeatCapacityRatioGamma(1.4);
	//configuration->setDistributionInitType(Iterative);


	parser.applyToSolverConfiguration(*configuration);


	boost::shared_ptr<ProblemDescription<2> > couetteProblem = cylinder;

    if(static_cast<bool>(parser.getArgument<int>("compressible"))!=true) {
        CFDSolver<2> solver(configuration, couetteProblem);
        solver.run();
    }
    else
    {
        CompressibleCFDSolver<2> solver(configuration, couetteProblem);
        solver.run();
    }

	pout << "NATriuM step-9 terminated." << endl;

	return 0;
}
