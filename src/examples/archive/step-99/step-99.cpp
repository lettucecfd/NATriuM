/*
 * step-99.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: kk
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"

#include "DensityFluctuation2D.h"

using namespace natrium;

// Main function
int main() {

	MPIGuard::getInstance();

	pout << "Starting NATriuM step-99 ..." << endl;

	size_t refinementLevel = 6;
	size_t orderOfFiniteElement = 2;
	bool isUnstructured = false;

	const double dqScaling = 1;
	const double viscosity = 0.003 ;
	const double startTime = 0.0;
	const double CFL = 0.4;
	// setup configuration

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<SolverConfiguration>();
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/step-99";
	configuration->setOutputDirectory(dirName.str());
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(1);
	configuration->setOutputTableInterval(100);
	configuration->setSimulationEndTime(5);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setCFL(CFL);
	configuration->setCommandLineVerbosity(7);

	//configuration->setAdvectionScheme(SEDG) ;
	//configuration->setCollisionScheme(BGK_WITH_TRANSFORMED_DISTRIBUTION_FUNCTIONS) ;
	//configuration->setStencil(Stencil_D2Q9) ;

	boost::shared_ptr<DensityFluctuation2D> DensityFluctuationFlow = boost::make_shared<DensityFluctuation2D>(viscosity, refinementLevel);
	boost::shared_ptr<ProblemDescription<2> > DensityFluctuationProblem = DensityFluctuationFlow;
	CFDSolver<2> solver(configuration, DensityFluctuationProblem);

	solver.run();

	pout << "step-99 terminated." << endl;

	return 0;
}
