/**
 * @file step-1.cpp
 * @short First tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */


#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-1..." << endl;

	size_t refinementLevel = 5;
	size_t orderOfFiniteElement = 3;

	shared_ptr<SolverConfiguration> configuration = make_shared<SolverConfiguration>();
	double deltaX = 1./(pow(2,refinementLevel)*(configuration->getOrderOfFiniteElement()-1));
	configuration->setOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setTimeStep(0.005*deltaX);
	configuration->setNumberOfTimeSteps(5000);
	configuration->setOutputDirectory( "../results/step-1-2");

	// set viscosity so that tau = 1
	double viscosity = 0.1*configuration->getTimeStep()/3.;

	shared_ptr<ProblemDescription<2> > taylorGreen = make_shared<TaylorGreenVortex2D>(viscosity, refinementLevel);

	CFDSolver<2> solver(configuration, taylorGreen);

	solver.run();

	cout << "NATriuM step-1 terminated." << endl;

	return 0;
}
