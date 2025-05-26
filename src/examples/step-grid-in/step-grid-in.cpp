/**
 * @date 31.03.2014
 * @author Andreas Kraemer, Dominik Wilde, Philipp Spelten, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

//! [Includes]
#include <stdlib.h>
#include <sstream>
#include <ctime>

#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D2Q25H.h"
#include "natrium/stencils/D2Q19V.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CommandLineParser.h"


#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "DiamondObstacle2D.h"
//! [Includes]

//! [Namespace]
using namespace natrium;
//! [Namespace]

//! [Main function]
int main(int argc, char** argv) {

    MPIGuard::getInstance();
    CommandLineParser parser(argc, argv);
    parser.setArgument<double>("Ma", "Mach number", 1.5);
    parser.setArgument<int>("Re", "Reynolds number", 10000);
    parser.setArgument<int>("ref-level", "Refinement level", 1);
    parser.setArgument<int>("aoa", "Angle of attack", 0);
    parser.setArgument<int>("server-end", "Maximum server time [s]", 70000);
    parser.setArgument<string>("foilname", "Folder of domain mesh", "nonUni_closed");

    try {
        parser.importOptions();
    } catch (HelpMessageStop&){
        return 0;
    }

    if (is_MPI_rank_0()) LOG(WELCOME) << "Starting NATriuM step-grid-in..." << endl;
    const int refLevel = parser.getArgument<int>("ref-level");

    // set Reynolds and Mach number
    const double Ma = parser.getArgument<double>("Ma")*sqrt(1.4);
    const double gamma = 1.4;
    // increase velocity to gain correct speed

    const double Re = parser.getArgument<int>("Re");

    // set Problem so that the right Re and Ma are achieved
    double U = 1;
    double reference_temperature = 1;
    double scaling = sqrt(3) * U / (Ma * sqrt(gamma*reference_temperature));
    const double viscosity = U / Re; // (because L = 1)
    int aoa = parser.getArgument<int>("aoa");

    // make problem and solver objects
	boost::shared_ptr<ProblemDescription<2>> obstacle_flow
        = boost::make_shared<DiamondObstacle2D>(U, viscosity, refLevel, aoa, parser.getArgument<string>("foilname"));
	//! [Problem]

	//! [Configuration]
	std::stringstream dirname;
    dirname << getenv("NATRIUM_HOME") << "/step-grid-in/Re" << Re << "-Ma" << Ma << "-reflevel" << refLevel
            << "-aoa" << aoa << "-time" << std::time(nullptr);
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<SolverConfiguration>();
	configuration->setOutputDirectory(dirname.str());
    configuration->setUserInteraction(false);
    configuration->setOutputCheckpointInterval(100000);
	configuration->setOutputSolutionInterval(10000);
    configuration->setStencilScaling(scaling);
	configuration->setNumberOfTimeSteps(200000);
    configuration->setHeatCapacityRatioGamma(gamma);
	//configuration->setTimeIntegrator(EXPONENTIAL);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
    configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
//	 configuration->setForcingScheme(NO_FORCING);
	configuration->setStencil(Stencil_D2Q19V);
	configuration->setMachNumber(Ma);
    configuration->setSupportPoints(GAUSS_LOBATTO_CHEBYSHEV_POINTS);
    configuration->setCollisionScheme(BGK_STANDARD);
    // configuration->setSimulationEndTime(30);
    configuration->setServerEndTime(parser.getArgument<int>("server-end"));

    parser.applyToSolverConfiguration(*configuration);

    CompressibleCFDSolver<2> solver(configuration, obstacle_flow);
    solver.run();

    if (is_MPI_rank_0()) LOG(WELCOME) << "NATriuM step-grid-in terminated." << endl;

	return 0;
}

//! [Solver]
