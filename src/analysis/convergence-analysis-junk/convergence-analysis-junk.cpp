/**
 * @file convergence-analysis-junk.cpp
 * @shortThe convergence of the NATriuM solver is analyzed by application to the Taylor-Green vortex in 2D (only periodic walls).
 * This script uses a diffusive scaling (= increasing Mach number). Thus, the numerical solution convergence against the
 * incompressible solution.
 * To analyze the results, move the table_order.txt and table_results.txt files to NATriuM/src/analysis/convergence_analysis_basic/
 * and execute the gnuplot scripts.
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/collision/BGKStandard.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

// Main function
int main() {

	MPIGuard::getInstance()

	pout << "Starting NATriuM convergence analysis (diffusive scaling)..."
			<< endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// Re = viscosity/(2*pi)
	const double viscosity = 1;
	//initial mach number = 0.05
	const double initialScaling = 100 / 5 * sqrt(3);
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = 2;
	//const double constant = -0.36;

	// prepare time table file
	std::stringstream filename;
	// the output is written to the standard output directory (e.g. NATriuM/results or similar)
	filename << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-junk/table_runtime.txt";
	std::ofstream timeFile(filename.str().c_str());
	timeFile
			<< "#refinement Level     dt        init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

	// prepare error table file
	std::stringstream filename2;
	filename2 << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-junk/table_order.txt";
	std::ofstream orderFile(filename2.str().c_str());
	orderFile << "# visc = " << viscosity << ";" << endl;
	orderFile
			<< "#  refinementlevel  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
			<< endl;

	for (size_t refinementLevel = 2; refinementLevel < 9; refinementLevel++) {
		pout << "refinement Level = " << refinementLevel << endl;
//		for (size_t orderOfFiniteElement = 2; orderOfFiniteElement < 7;
//				orderOfFiniteElement++) {
//			pout << "FE order = " << orderOfFiniteElement << endl;
		// the scaling has to be orders of magnitude greater than the boundary velocity

		// calculate distance between quadrature nodes
		double dx = 2 * 3.1415926
				/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));

		// chose scaling so that the ratio between log xi_0^2 and log tau is constant
		double scaling = initialScaling / sqrt(dx);	//pow(3*viscosity*sqrt(2)/dx,constant/(constant+2.0));

		// chose dt so that courant (advection) = 1 for the diagonal directions
		// double dt = dx / (scaling * sqrt(2));
		//chose CFL = 0.4
		double dt = 0.4 * dx / scaling;

		pout << "dt = " << dt;

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/convergence-analysis-junk/"
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
			pout << "Timestep too big." << endl;
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
		const double tau = BGK::calculateRelaxationParameter(
				viscosity, dt, *solver.getStencil());
		pout << "... scaling = " << scaling << " ... tau = " << tau << " ...";
		//pout << "constant = " << log2(scaling*scaling) / log2(tau) << " ...";
		time1 = clock() - timestart;

		try {
			solver.run();
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;
			pout << " OK ... Init: " << time1 << " sec; Run: " << time2
					<< " sec." << endl;
			// put out runtime
			timeFile << refinementLevel << "         " << dt << "      "
					<< time1 << "     " << time2 << "        "
					<< time2 / configuration->getNumberOfTimeSteps() << endl;
			// put out final errors
			solver.getErrorStats()->update();
			orderFile << refinementLevel << " " << solver.getIteration() << " "
					<< solver.getTime() << " "
					<< solver.getErrorStats()->getMaxUAnalytic() << " "
					<< solver.getErrorStats()->getMaxVelocityError() << " "
					<< solver.getErrorStats()->getMaxDensityError() << " "
					<< solver.getErrorStats()->getL2VelocityError() << " "
					<< solver.getErrorStats()->getL2DensityError() << endl;
		} catch (std::exception& e) {
			pout << " Error" << endl;
		}
		//}
		//}
//	}
	}
	pout << "Convergence analysis (diffusive scaling) terminated." << endl;

	return 0;
}
