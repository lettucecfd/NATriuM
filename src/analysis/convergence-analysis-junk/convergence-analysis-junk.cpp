/**
 * @file convergence-analysis-basic.cpp
 * @short Taylor-Green vortex in 2D (only periodic walls)
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "../../examples/step-1/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM convergence analysis (basic)..." << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	const double viscosity = 1;
	//initial mach number
	const double initialMa = 0.05;
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = 2;
	const double constant = -0.36;

	for (size_t refinementLevel = 2; refinementLevel < 10; refinementLevel++) {
		cout << "refinement Level = " << refinementLevel << endl;
//		for (size_t orderOfFiniteElement = 2; orderOfFiniteElement < 7;
//				orderOfFiniteElement++) {
//			cout << "FE order = " << orderOfFiniteElement << endl;
		// the scaling has to be orders of magnitude greater than the boundary velocity

		// calculate distance between quadrature nodes
		double dx = 2 * 3.1415926
				/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));

		// chose scaling so that the ratio between log xi_0^2 and log tau is constant
		double scaling = pow(3*viscosity*sqrt(2)/dx,constant/(constant+2.0));

		// chose dt so that courant (advection) = 1 for the diagonal directions
		double dt = dx / (scaling * sqrt(2));


		cout << "dt = " << dt << " ..." ;

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << "../results/convergence-analysis-junk/"
				<< orderOfFiniteElement << "_" << refinementLevel;
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(10);
		configuration->setOutputCheckpointInterval(1000);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(0);
		configuration->setTimeStepSize(dt);
		if (dt > 0.1) {
			cout << "Timestep too big." << endl;
			continue;

		}
		configuration->setNumberOfTimeSteps(1.0 / dt);

		// make problem and solver objects
		shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel);
		shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
		timestart = clock();
		BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
		// get tau to test if the "constant value" is really constant
		const double tau = BGKTransformed::calculateRelaxationParameter(
						viscosity, dt, solver.getBoltzmannModel());
		cout << "constant = " << log2(scaling*scaling) / log2(tau) << " ...";
		time1 = clock() - timestart;

		try {
			solver.run();
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;
			cout << " OK ... Init: " << time1 << " sec; Run: " << time2
					<< " sec." << endl;
		} catch (std::exception& e) {
			cout << " Error" << endl;
		}
		//}
		//}
//	}
	}
	cout << "Convergence analysis (basic) terminated." << endl;

	return 0;
}