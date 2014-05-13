/**
 * @file step-2.cpp
 * @short Second tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "deal.II/numerics/data_out.h"

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "CouetteFlow2D.h"

//#define PRINT_SYSTEM_VECTOR

using namespace natrium;

void getAnalyticSolution(double time, distributed_vector& analyticSolution1,
		distributed_vector& analyticSolution2,
		const vector<dealii::Point<2> >& supportPoints,
		const CouetteFlow2D& couette) {
	assert(analyticSolution1.size() == supportPoints.size());
	assert(analyticSolution2.size() == supportPoints.size());
	assert(supportPoints.size() > 0);

	for (size_t i = 0; i < supportPoints.size(); i++) {
		analyticSolution1(i) = couette.analyticVelocity1(supportPoints.at(i),
				time);
		analyticSolution2(i) = couette.analyticVelocity2(supportPoints.at(i),
				time);
	}
}

// Main function
int main() {

	cout << "Starting NATriuM step-2..." << endl;

	// set parameters, set up configuration object
	size_t refinementLevel = 4;
	size_t orderOfFiniteElement = 5;
	const double dqScaling = 1; //2 * sqrt(3);

	// chose U (the velocity of the top wall) so that Ma = 0.05
	const double U = 5. / 100. * sqrt(3) * dqScaling;
	// chose viscosity so that Re = 2000
	const double viscosity = U / 2000.;

	cout << "Mach number: " << U / (sqrt(3) * dqScaling) << endl;
	// configure solver
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory("../results/step-2");
	configuration->setRestart(false);
	configuration->setOutputFlags(
			configuration->getOutputFlags() | out_Checkpoints);
	configuration->setOutputCheckpointEvery(1000);
	configuration->setOutputVectorFieldsEvery(100000000);
	configuration->setNumberOfTimeSteps(100000000);
	configuration->setOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setDQScaling(dqScaling);
	configuration->setTimeStep(0.00025);
	//configuration->setDistributionInitType(Iterative);

	shared_ptr<CouetteFlow2D> couetteFlow = make_shared<CouetteFlow2D>(
			viscosity, U, refinementLevel);
	shared_ptr<ProblemDescription<2> > couetteProblem = couetteFlow;
	CFDSolver<2> solver(configuration, couetteProblem);
	cout << "Number of DoFs: " << 9*solver.getNumberOfDoFs() << endl;

#ifdef PRINT_SYSTEM_VECTOR
	// put out system vector
	std::stringstream s1;
	s1 << configuration->getOutputDirectory().c_str() << "/system_vector.txt";
	std::ofstream vectorOut(s1.str().c_str());
	solver.getAdvectionOperator()->getSystemVector().print(vectorOut);
	vectorOut << endl << endl << endl;
	solver.getAdvectionOperator()->getSystemMatrix().print_formatted(vectorOut);
	vectorOut.close();
#endif

	// File for max errors
	std::stringstream s;
	s << configuration->getOutputDirectory().c_str() << "/max_error.txt";
	std::fstream* maxNormOut;
	if (solver.getIterationStart() > 0) {
		maxNormOut = new std::fstream(s.str().c_str(),
				std::fstream::out | std::fstream::app);
	} else {
		maxNormOut = new std::fstream(s.str().c_str(), std::fstream::out);
		(*maxNormOut)
				<< "# t           u_numeric(at) y=L/2       max|u_numeric- u_analytic|    max |rho-rho_0|"
				<< endl;
	}

	// THE LOOP
	// The easiest thing at this point would be to simply type
	//solver.run();
	// but we want to do more than just the standard CFD simulation:
	// We want to calculate the error between the numerical and the analytic solution of the flow problem.
	// We need the support points to relate degrees of freedom to space coordinates.

	vector<dealii::Point<2> > supportPoints(
			solver.getAdvectionOperator()->getDoFHandler()->n_dofs());
	dealii::DoFTools::map_dofs_to_support_points(
			solver.getAdvectionOperator()->getMapping(),
			*solver.getAdvectionOperator()->getDoFHandler(), supportPoints);
	//find out middle of domain
	size_t middleDof = 0;
	double distToMiddle = 1000000.0;
	for (size_t i = 0; i < solver.getNumberOfDoFs(); i++) {
		double dist = pow(0.5 - supportPoints.at(i)(0), 2)
				+ pow(0.5 - supportPoints.at(i)(1), 2);
		if (dist < distToMiddle) {
			dist = distToMiddle;
			middleDof = i;
		}
	}
	distributed_vector analyticSolution1(solver.getNumberOfDoFs());
	distributed_vector analyticSolution2(solver.getNumberOfDoFs());
	double t = 0;
	for (size_t i = solver.getIterationStart(); t < 40; i++) {
		t = i * configuration->getTimeStep();
		if (i % 1000 == 0) {
			cout << "Iteration " << i << ",  t = " << t << endl;
		}
		// Stream and collide
		solver.output(i);
		solver.stream();
		solver.collide();

		if (i % 1 == 0) {
			if (t > 0.001) {
				// put out max velocity norm for numerical and analytic solution
				getAnalyticSolution(configuration->getTimeStep() * (i + 1),
						analyticSolution1, analyticSolution2, supportPoints,
						*couetteFlow);

				/// vtu ///
				std::stringstream str;
				str << configuration->getOutputDirectory().c_str() << "/t_" << i
						<< ".vtu";
				std::ofstream outfile(str.str().c_str());
				dealii::DataOut<2> data_out;
				data_out.attach_dof_handler(
						*solver.getAdvectionOperator()->getDoFHandler());
				data_out.add_data_vector(solver.getDensity(), "rho");
				data_out.add_data_vector(solver.getVelocity().at(0), "vx");
				data_out.add_data_vector(solver.getVelocity().at(1), "vy");
				// calculate analytic solution
				getAnalyticSolution(configuration->getTimeStep() * i,
						analyticSolution1, analyticSolution2, supportPoints,
						*couetteFlow);
				data_out.add_data_vector(analyticSolution1, "vx_analytic");
				data_out.add_data_vector(analyticSolution2, "vy_analytic");
				data_out.build_patches();
				data_out.write_vtu(outfile);

				// max error file //
				// calculate error
				analyticSolution1 *= -1.0;
				analyticSolution1.add(solver.getVelocity().at(0));
				analyticSolution2 *= -1.0;
				analyticSolution2.add(solver.getVelocity().at(1));
				double xVelocityInTheMiddle = solver.getVelocity().at(0)(
						middleDof);
				(*maxNormOut) << i * configuration->getTimeStep() << "  "
						<< xVelocityInTheMiddle << "  "
						<< analyticSolution1.linfty_norm() << " "
						<< solver.getMaxDensityDeviationFrom(1) << endl;
			}
		}
	}
	delete maxNormOut;

	cout << "NATriuM step-2 terminated." << endl;

	return 0;
}
