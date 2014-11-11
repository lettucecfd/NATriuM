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

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "DensityFluctuation2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-99 ..." << endl;

	size_t refinementLevel = 6;
	size_t orderOfFiniteElement = 2;
	bool isUnstructured = false;

	const double dqScaling = 1;
	const double viscosity = 0.003 ;
	const double startTime = 0.0;
	const double timeStepSize = 0.01;
	// setup configuration

	shared_ptr<SolverConfiguration> configuration = make_shared<SolverConfiguration>();
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/step-99";
	configuration->setOutputDirectory(dirName.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(1);
	configuration->setOutputTableInterval(100);
	configuration->setNumberOfTimeSteps(5./timeStepSize);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(dqScaling);
	configuration->setTimeStepSize(timeStepSize);
	configuration->setCommandLineVerbosity(7);

	//configuration->setAdvectionScheme(SEDG) ;
	//configuration->setCollisionScheme(BGK_WITH_TRANSFORMED_DISTRIBUTION_FUNCTIONS) ;
	//configuration->setStencil(Stencil_D2Q9) ;

	shared_ptr<DensityFluctuation2D> DensityFluctuationFlow = make_shared<DensityFluctuation2D>(viscosity, refinementLevel);
	shared_ptr<ProblemDescription<2> > DensityFluctuationProblem = DensityFluctuationFlow;
	CFDSolver<2> solver(configuration, DensityFluctuationProblem);

	solver.run();

	cout << "step-99 terminated." << endl;

	return 0;
}
