/*
 * step-6.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: bajat
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>




#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/Stencil.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/CommandLineParser.h"

#include "natrium/benchmarks/TaylorGreenVortex3D.h"


using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-17 ..." << endl;

	/////////////////////////////////////////////////
	// read from command line
	//////////////////////////////////////////////////

	CommandLineParser parser(argc, argv);
	parser.setArgument<int>("Re", "Reynolds number 1/nu", 1600);
    parser.setArgument<double>("Ma", "Mach number", 0.1);
    parser.setArgument<int>("densityNumerator", "Is the numerator in the density dependent on the Mach number?", 0);

    parser.setPositionalArgument<int>("ref-level",
			"Refinement level of the computation grid.");
	parser.setArgument<int>("grid-repetitions",
			"Number of grid cells along each axis before global refinement; "
			"to produce grids with refinements that are not powers of two.", 1);
	try {
		parser.importOptions();
	} catch (HelpMessageStop&) {
		return 0;
	}
	double Re = parser.getArgument<int>("Re");
	double refinement_level = parser.getArgument<int>("ref-level");
	double repetitions = parser.getArgument<int>("grid-repetitions");
    const bool isDensityNumerator = static_cast<bool>(parser.getArgument<int>("densityNumerator"));
    const double densityNumerator = isDensityNumerator ? parser.getArgument<double>("Ma") *
                                                                                     parser.getArgument<double>("Ma") *
                                                                                     1.4 : 1.0;
	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// im Paper von Gassner und Beck ist U = 1/2pi definiert !!!!!
	// Aber ihre zeitangaben beziehen sich auf U = 1, wie bei Brachet (1991)
	// hier simulieren wir jetzt  U = 1 (ist im TGV3D modul sowieso nur so definiert)
	// und mit der Reynoldsähnlichkeit ist das kein Problem, da wir gegenüber Gassner
	// und Beck ja auch die Viskosität verändern
	// Nach einem Blick in van Rees et.al. (2011) und einer erfolgreichen Simulation in Palabos: 
	// (Hier war tau = 3*nu_LB + 0.5, nu_LB = U_lattice * (N/2pi) / Re, dt = 2pi/N* U_lattice
	// Re = 1/nu,L=2pi, U = 1 und D = [0,2pi*L]^3
	const double U = 1;
	//const double L = 2 * M_PI;
	const double viscosity = 1.0 / Re;
    const double Ma = parser.getArgument<double>("Ma")*sqrt(1.4);
	const double cs = U / Ma;

	// chose scaling so that the right Ma-number is achieved
	const double scaling = sqrt(3) * cs;
	const bool init_rho_analytically = true;


	boost::shared_ptr<ProblemDescription<3> > taylorGreen = boost::make_shared<
			TaylorGreenVortex3D>(viscosity, refinement_level, cs,
			init_rho_analytically, repetitions, densityNumerator);

	// setup configuration
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setUserInteraction(false);
	configuration->setOutputCheckpointInterval(1e9);
	configuration->setOutputSolutionInterval(10000);
	configuration->setSimulationEndTime(10.0);
	configuration->setOutputGlobalTurbulenceStatistics(false);
    configuration->setOutputCompressibleTurbulenceStatistics(true);
	configuration->setStencilScaling(scaling);
	configuration->setStencil(Stencil_D3Q45);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
	configuration->setHeatCapacityRatioGamma(1.4);
	configuration->setPrandtlNumber(0.71);

	parser.applyToSolverConfiguration(*configuration);

	// standard output dir
	if (not parser.hasArgument("output-dir")){
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/step-17-compressible/Re" << Re << "-ref"
				<< refinement_level << "-p"
				<< configuration->getSedgOrderOfFiniteElement() << "-coll"
				<< static_cast<int>(configuration->getCollisionScheme()) << "-sl"
				<< static_cast<int>(configuration->getAdvectionScheme());
		if (configuration->getAdvectionScheme() != SEMI_LAGRANGIAN)
			dirName << "-int"
					<< static_cast<int>(configuration->getTimeIntegrator()) << "_"
					<< static_cast<int>(configuration->getDealIntegrator());
		dirName << "-CFL" << configuration->getCFL();
		dirName << "-sten" << static_cast<int>(configuration->getStencil());
        dirName << "-num" << static_cast<double>(densityNumerator);
		if (Ma!=0.1)
		    dirName << "-Ma" << Ma;
		if (configuration->isFiltering())
			dirName << "-filt"
					<< static_cast<int>(configuration->getFilteringScheme())
					<< "by_max_degree";
		if (configuration->getRegularizationScheme() != NO_REGULARIZATION)
			dirName << "-reg"
					<< static_cast<int>(configuration->getRegularizationScheme());
        if (configuration->getEquilibriumScheme()!= BGK_EQUILIBRIUM)
            dirName << "-equili"
                    << static_cast<int>(configuration->getEquilibriumScheme());
		if (configuration->getCollisionScheme() == MRT_STANDARD) {
			dirName << "-mrt" << static_cast<int>(configuration->getMRTBasis());
		}
		if (configuration->getCollisionScheme() == MRT_STANDARD) {
			dirName << "-relax"
					<< static_cast<int>(configuration->getMRTRelaxationTimes());
		}
		if (repetitions != 1) {
					dirName << "-rep" << repetitions;
				}
		configuration->setOutputDirectory(dirName.str());
	}
	/////////////////////////////////////////////////
	// run solver
	//////////////////////////////////////////////////

	CompressibleCFDSolver<3> solver(configuration, taylorGreen);

	const size_t table_output_lines_per_10s = 300;
	configuration->setOutputTableInterval(
			1 + 10.0 / solver.getTimeStepSize() / table_output_lines_per_10s);

	solver.run();

	pout << "step-17 terminated." << endl;

	return 0;

}
