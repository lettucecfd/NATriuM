/**
 * @file step-3.cpp
 * @short Third Tutorial:  Lid-Driven cavitys in 2D
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "LidDrivenCavity2D.h"

using namespace natrium;

// analytic solution (where should the bump be after a certain time)
void getAnalyticSolution(double time, distributed_vector& analyticSolution1,
		distributed_vector& analyticSolution2,
		const vector<dealii::Point<2> >& supportPoints,
		const LidDrivenCavity2D& ldCavity) {
	assert(analyticSolution1.size() == supportPoints.size());
	assert(analyticSolution2.size() == supportPoints.size());
	assert(supportPoints.size() > 0);

	for (size_t i = 0; i < supportPoints.size(); i++) {
		analyticSolution1(i) = ldCavity.analyticVelocity1(supportPoints.at(i),
				time);
		analyticSolution2(i) = ldCavity.analyticVelocity2(supportPoints.at(i),
				time);
	}
}

// Main function
int main() {

	cout << "Starting NATriuM step-3..." << endl;

	// set Reynolds and Mach number
	const double Re = 1000;
	const double Ma = 0.1;

	// set spatial discretization
	size_t refinementLevel = 4;
	size_t orderOfFiniteElement = 5;

	// set temporal discretization
	double CFL = 0.1;
	double timestep = CFL * 1./pow(2,4) /sqrt(2);

	// set Problem so that the right Re and Ma are achieved
	const double U = 1 / sqrt(3) * Ma;
	const double viscosity = U / Re; // (because L = 1)

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-3";
	//dirname << getenv("NATRIUM_HOME") << "/step-3-implicit";
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(1000);
	//configuration->setTimeIntegrator(THETA_METHOD);
	configuration->setOutputSolutionInterval(100);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(1);
	configuration->setTimeStepSize(timestep);
	configuration->setNumberOfTimeSteps(200000);
	configuration->setInitializationScheme(EQUILIBRIUM);

	// make problem and solver objects
	shared_ptr<LidDrivenCavity2D> lidDrivenCavity = make_shared<
			LidDrivenCavity2D>(U, viscosity, refinementLevel);
	shared_ptr<ProblemDescription<2> > ldCavityProblem = lidDrivenCavity;
	CFDSolver<2> solver(configuration, ldCavityProblem);

	solver.run();

	cout << "NATriuM step-3 terminated." << endl;

	return 0;
}
