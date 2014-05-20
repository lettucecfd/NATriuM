/**
 * @file step-3.cpp
 * @short Third Tutorial:  Lid-Driven cavitys in 2D
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

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

	// set parameters, set up configuration object
	size_t refinementLevel = 5;
	size_t orderOfFiniteElement = 2;
	double viscosity = 0.0001/(sqrt(3));

	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	double deltaX = 1.
			/ (pow(2, refinementLevel)
					* (configuration->getSedgOrderOfFiniteElement() - 1));
	configuration->setOutputDirectory("../results/step-3");
	configuration->setRestartAtLastCheckpoint(true);
	configuration->setOutputCheckpointInterval(1000);
	configuration->setOutputSolutionInterval(100);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(50);
	double tScaling =  0.4 / (sqrt(2) * configuration->getStencilScaling());
	configuration->setTimeStepSize(tScaling * deltaX);
	configuration->setTimeStepSize( 1.2e-05);
	configuration->setNumberOfTimeSteps(100000);
	configuration->setInitializationScheme(ITERATIVE);

	// make problem and solver objects
	shared_ptr<LidDrivenCavity2D> lidDrivenCavity = make_shared<LidDrivenCavity2D>(
			viscosity, refinementLevel);
	shared_ptr<ProblemDescription<2> > ldCavityProblem = lidDrivenCavity;
	CFDSolver<2> solver(configuration, ldCavityProblem);

	// File for max norms
	/*	std::stringstream s;
	 s << configuration->getOutputDirectory().c_str() << "/max_norm.txt";
	 std::fstream* maxNormOut;
	 if (solver.getIterationStart() > 0) {
	 maxNormOut = new std::fstream(s.str().c_str(),
	 std::fstream::out | std::fstream::app);
	 } else {
	 maxNormOut = new std::fstream(s.str().c_str(), std::fstream::out);
	 (*maxNormOut)
	 << "# t            max|u_numeric|            max|u_analytic|     max |rho-rho_0|"
	 << endl;
	 }
	 */

	solver.run();
	// THE LOOP
	// The easiest thing at this point would be to simply type
	//solver.run();
	// but we want to do more than just the standard CFD simulation:
	// We want to calculate the error between the numerical and the analytic solution of the flow problem.
	// We need the support points to relate degrees of freedom to space coordinates.
	/*	vector<dealii::Point<2> > supportPoints(
	 solver.getAdvectionOperator()->getDoFHandler()->n_dofs());
	 dealii::DoFTools::map_dofs_to_support_points(
	 solver.getAdvectionOperator()->getMapping(),
	 *solver.getAdvectionOperator()->getDoFHandler(), supportPoints);
	 distributed_vector analyticSolution1(solver.getNumberOfDoFs());
	 distributed_vector analyticSolution2(solver.getNumberOfDoFs());
	 size_t N = configuration->getNumberOfTimeSteps();
	 for (size_t i = solver.getIterationStart(); i < N; i++) {
	 if (i % 100 == 0) {
	 cout << "Iteration " << i << endl;
	 }
	 // Stream and collide
	 solver.output(i);
	 solver.stream();
	 solver.collide();

	 //custom output
	 std::stringstream str;
	 str << configuration->getOutputDirectory().c_str() << "/t_" << i
	 << ".vtu";
	 std::string filename = str.str();
	 std::ofstream vtu_output(filename.c_str());
	 dealii::DataOut<2> data_out;
	 data_out.attach_dof_handler(
	 *solver.getAdvectionOperator()->getDoFHandler());
	 data_out.add_data_vector(solver.getDensity(), "rho");
	 data_out.add_data_vector(solver.getVelocity().at(0), "v_1");
	 data_out.add_data_vector(solver.getVelocity().at(1), "v_2");
	 // calculate analytic solution
	 getAnalyticSolution(configuration->getTimeStep() * i, analyticSolution1,
	 analyticSolution2, supportPoints, *tgVortex);
	 data_out.add_data_vector(analyticSolution1, "v_1_analytic");
	 data_out.add_data_vector(analyticSolution2, "v_2_analytic");
	 data_out.build_patches();
	 data_out.write_vtu(vtu_output);
	 // put out max velocity norm for numerical and analytic solution
	 (*maxNormOut) << i * configuration->getTimeStep() << "  "
	 << solver.getMaxVelocityNorm() << "  "
	 << exp(-2 * viscosity * i * configuration->getTimeStep()) << " "
	 << solver.getMaxDensityDeviationFrom(1) << endl;
	 }
	 delete maxNormOut;

	 */

	cout << "NATriuM step-3 terminated." << endl;

	return 0;
}
