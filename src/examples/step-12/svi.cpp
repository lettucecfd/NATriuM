/**
 * @file static_droplet.cpp
 * @short Static droplet simulation
 * @date 16.11.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"
#include "deal.II/base/utilities.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/ProblemDescription.h"


#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"

#include "ShockVortexInteraction.h"

#include "natrium/utilities/Info.h"
#include "natrium/utilities/CommandLineParser.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// ========================================================================
	// READ COMMAND LINE PARAMETERS
	// ========================================================================

	CommandLineParser parser(argc, argv);
    parser.addDocumentationString("svi",
			"2D Kelvin-Helmholtz instability at Re=30'000, as in Boesch (2014)");
	parser.setPositionalArgument<int>("ref-level",
			"refinement of the computational grid");
	parser.setArgument<double>("tx",
			"transformation of the grid in x-direction (<1)", 0);
	parser.setArgument<double>("ty",
			"transformation or the grid in y-direction (<1)", 0);
	parser.setArgument<int>("filter", "apply filtering", 0);
	parser.setArgument<int>("filter-s", "parameter as filter", 32);
    parser.setArgument<int>("vmult", "apply vMultLimiter", 0);
    parser.setFlag("minion-brown",
            "sets the problem up as the 'thin' shear-layer in the original work by Minion and Brown");
	try {
		parser.importOptions();
	} catch (HelpMessageStop&) {
		return 0;
	}

	// ========================================================================
	// MAKE FLOW PROBLEM
	// ========================================================================
    //Speed of the vortex
	double Ma_v = 0.5;
	double perturbation = 0.05;
	double kappa = 80;
    double u0 = 0.0;

    double t_max=1.0;

    double scaling = 1.0;
    double viscosity = 0.001725164;
    boost::shared_ptr<ProblemDescription<2> > svi = boost::make_shared<
            ShockVortexInteraction>(viscosity, parser.getArgument<int>("ref-level"), u0,
			kappa, Ma_v, perturbation, parser.getArgument<double>("tx"),
			parser.getArgument<double>("ty"));

	// ========================================================================
	// CONFIGURE SOLVER
	// ========================================================================
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setSwitchOutputOff(false);
	configuration->setUserInteraction(false);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setOutputTableInterval(10);	//10
	configuration->setOutputSolutionInterval(10); //10
	configuration->setOutputCheckpointInterval(1e9);
	configuration->setConvergenceThreshold(1e-10);
	configuration->setStencilScaling(scaling);
	configuration->setCFL(1);
	configuration->setSedgOrderOfFiniteElement(2);
	configuration->setSimulationEndTime(t_max);
	configuration->setOutputGlobalTurbulenceStatistics(true);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	configuration->setStencil(Stencil_D2Q25H);
	configuration->setExponentialFilterAlpha(36);
	configuration->setExponentialFilterNc(4);
	configuration->setPrandtlNumber(0.75);
	configuration->setHeatCapacityRatioGamma(1.4);
	configuration->setInitializationScheme(EQUILIBRIUM);
	configuration->setIterativeInitializationNumberOfIterations(300);
    configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);


    configuration->setVmultLimiter(bool(parser.getArgument<int>("vmult")));
    pout << "VMultLimiter is " << configuration->isVmultLimiter() << endl;

	parser.applyToSolverConfiguration(*configuration);


	std::stringstream dirname;
    dirname << getenv("NATRIUM_HOME") << "/ShockVortexInteraction";

	dirname << "/N" << parser.getArgument<int>("ref-level") << "-p"
			<< configuration->getSedgOrderOfFiniteElement() << "-sl"
			<< static_cast<int>(configuration->getAdvectionScheme()) << "-coll"
			<< static_cast<int>(configuration->getCollisionScheme()) << "-int"
			<< static_cast<int>(configuration->getTimeIntegrator()) << "_"
			<< static_cast<int>(configuration->getDealIntegrator()) << "-CFL"
			<< configuration->getCFL() << "-reg" << static_cast<int>(configuration->getRegularizationScheme())<< "-scaling"
            << configuration->getStencilScaling() << "-suppP"
            << configuration->getSupportPoints() << "-Pr"
            << configuration->getPrandtlNumber() << "-sten"
            << configuration->getStencil() << "-Ma_v"
            << Ma_v ;
	if (parser.getArgument<int>("filter") != 0) {
		dirname << "-filter" << parser.getArgument<int>("filter") << "-filt_s"
				<< parser.getArgument<int>("filter-s");
	}
	if ((parser.getArgument<double>("tx") != 0)
			or (parser.getArgument<double>("ty") != 0)) {
		dirname << "-tx" << parser.getArgument<double>("tx") << "-ty"
				<< parser.getArgument<double>("ty");
	}
	if (configuration->getRegularizationScheme() != NO_REGULARIZATION){
		dirname << "-reg" << static_cast<int>(configuration->getRegularizationScheme());
	}
	configuration->setOutputDirectory(dirname.str());

	pout << "Simulation end time will be t_max = " << t_max << endl;
	// ========================================================================
	// RUN SOLVER
	// ========================================================================

    natrium::CompressibleCFDSolver<2> solver(configuration, svi);

	solver.run();

	// ========================================================================
	// FINAL OUTPUT
	// ========================================================================
	pout << "Simulation successful." << endl;
	return 0;
}
