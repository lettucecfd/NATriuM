/**
 * @file step-1.cpp
 * @short First tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-1..." << endl;

	// set parameters, set up configuration object
	//size_t refinementLevel = 5;
	//size_t orderOfFiniteElement = 2;
	double viscosity = 1;

	cout << "Step-1 " << endl;
	for (size_t scaling = 256; scaling >= 1; scaling /= 4) {
		cout << "Stencil scaling = " << scaling << endl;
		for (size_t refinementLevel = 2; refinementLevel < 5;
							refinementLevel++) {
						cout << "refinement Level = " << refinementLevel << endl;
for (size_t orderOfFiniteElement = 6; orderOfFiniteElement < 10;
				orderOfFiniteElement++) {
			cout << "FE order = " << orderOfFiniteElement << endl;
							for (double dt = 0.001; dt >= 0.0001; dt /= 2) {
					cout << "dt = " << dt << " ...";

					// time measurement variables
					double time1, time2, timestart;

					std::stringstream dirName;
					dirName << "../results/step-1/" << orderOfFiniteElement
							<< "_" << refinementLevel << "_" << scaling << "_"
							<< dt;

					shared_ptr<SolverConfiguration> configuration = make_shared<
							SolverConfiguration>();
					configuration->setOutputDirectory(dirName.str().c_str());
					configuration->setRestartAtLastCheckpoint(false);
					configuration->setUserInteraction(false);
					configuration->setOutputTableInterval(100);
					configuration->setOutputCheckpointInterval(1000);
					configuration->setSedgOrderOfFiniteElement(
							orderOfFiniteElement);
					configuration->setStencilScaling(scaling);
					configuration->setCommandLineVerbosity(0);
					//double tScaling = std::min(0.1,
					//		1. / (2 * configuration->getStencilScaling()));
					//double deltaX =
					//		1.
					//				/ (pow(2, refinementLevel)
					//						* (configuration->getSedgOrderOfFiniteElement()
					//								- 1));
					//configuration->setTimeStepSize(tScaling * deltaX);
					configuration->setTimeStepSize(dt);
					configuration->setNumberOfTimeSteps(0.5 / dt);

					//configuration->setDistributionInitType(Iterative);

					// make problem and solver objects
					shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
							TaylorGreenVortex2D>(viscosity, refinementLevel);
					shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
					timestart = clock();
					BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
					time1 = clock() -timestart;

					try {
						solver.run();
						time2 = clock() - time1 - timestart;
						time1 /= CLOCKS_PER_SEC;
						time2 /= CLOCKS_PER_SEC;
						cout << " OK ... Init: " << time1 << " sec; Run: " << time2 << " sec." << endl;
					} catch (std::exception& e) {
						cout << " Error" << endl;
					}
				}
			}
		}
	}
	cout << "NATriuM step-1 terminated." << endl;

	return 0;
}
