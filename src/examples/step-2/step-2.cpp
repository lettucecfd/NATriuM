/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */


#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "CouetteFlow2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-2..." << endl;

	// set parameters, set up configuration object
	size_t refinementLevel = 4;
	size_t orderOfFiniteElement = 2;
	double viscosity = 0.01;

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	double deltaX = 1.
			/ (pow(2, refinementLevel)
					* (configuration->getOrderOfFiniteElement() - 1));
	configuration->setOutputDirectory("../results/step-2");
	configuration->setRestart(false);
	configuration->setOutputFlags(
			configuration->getOutputFlags() | out_Checkpoints);
	configuration->setOutputCheckpointEvery(500);
	configuration->setOutputVectorFieldsEvery(10);
	configuration->setOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setDQScaling(2*sqrt(3));
	double tScaling = std::min(0.1, 1. / (2 * configuration->getDQScaling()));
	configuration->setTimeStep(0.001);
	configuration->setNumberOfTimeSteps(50000);
	//configuration->setDistributionInitType(Iterative);

	double topPlateVelocity = 0.1;
	shared_ptr<ProblemDescription<2> > couetteFlow = make_shared<CouetteFlow2D>(viscosity, topPlateVelocity, refinementLevel);
	CFDSolver<2> solver(configuration, couetteFlow);

	solver.run();

	cout << "NATriuM step-2 terminated." << endl;

	return 0;
}
