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

#include "CouetteFlow2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-1..." << endl;

	double relaxationParameter = 0.7;
	double topPlateVelocity = 0.05;

	shared_ptr<ProblemDescription<2> > couetteFlow = make_shared<CouetteFlow2D>(relaxationParameter, topPlateVelocity);
	shared_ptr<SolverConfiguration> configuration = make_shared<SolverConfiguration>();
	CFDSolver<2> solver(configuration, couetteFlow);

	solver.run();

	cout << "NATriuM step-1 terminated." << endl;

	return 0;
}
