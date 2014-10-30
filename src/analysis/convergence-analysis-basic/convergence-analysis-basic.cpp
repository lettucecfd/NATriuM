/**
 * @file convergence-analysis-basic.cpp
 * @short The convergence of the NATriuM solver is analyzed by application to the Taylor-Green vortex in 2D (only periodic walls).
 * This script uses a linear scaling (= constant Mach number). Thus, there is a general compressibily error, which destroys the
 * convergence for the finer meshes. To analyze the results, move the table_order.txt and table_results.txt files to NATriuM/src/analysis/convergence_analysis_basic/
 * and execute the gnuplot scripts.
 * @date 05.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "solver/BenchmarkCFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/Benchmark.h"

#include "utilities/BasicNames.h"

#include "../../examples/step-1/TaylorGreenVortex2D.h"

using namespace natrium;

// if this define statement is enabled: only the initialization time is regarded
//#define MEASURE_ONLY_INIT_TIME

// Main function
int main() {

	cout << "Starting NATriuM convergence analysis (linear scaling)..." << endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	const double viscosity = 1;
	// C-E-approach: constant stencil scaling
	// specify Mach number
	const double Ma = 0.05;
	// zunaechst: fixed order of FE
	const double orderOfFiniteElement = 2;

	// chose scaling so that the right Mach number is achieved
	double scaling = sqrt(3) * 1 / Ma;

	// prepare time table file
	// the output is written to the standard output directory (e.g. NATriuM/results or similar)
	std::stringstream filename;
	filename << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-basic/table_runtime.txt";
	std::ofstream timeFile(filename.str().c_str());
	timeFile
			<< "#refinement Level     dt        init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

	// prepare error table file
	std::stringstream filename2;
	filename2 << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-basic/table_order.txt";
	std::ofstream orderFile(filename2.str().c_str());
	orderFile << "# visc = " << viscosity << "; Ma = " << Ma << endl;
	orderFile
			<< "#  refinementlevel  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
			<< endl;

	for (size_t refinementLevel = 2; refinementLevel < 9; refinementLevel++) {
		cout << "refinement Level = " << refinementLevel << endl;

		double dx = 2 * 3.1415926
				/ (pow(2, refinementLevel) * (orderOfFiniteElement - 1));
		// chose dt so that courant (advection) = 1 for the diagonal directions
		//double dt = dx / (scaling * sqrt(2));
		//chose CFL = 0.4
		double dt = 0.4 * dx / scaling;

		cout << "dt = " << dt << " ...";

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/convergence-analysis-basic/"
				<< orderOfFiniteElement << "_" << refinementLevel;
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		//configuration->setOutputTableInterval(10);
		//configuration->setOutputCheckpointInterval(1000);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(0);
		configuration->setTimeStepSize(dt);
		if (dt > 0.1) {
			cout << "Timestep too big." << endl;
			continue;

		}
		configuration->setNumberOfTimeSteps(1.0 / dt);

#ifdef MEASURE_ONLY_INIT_TIME
		configuration->setNumberOfTimeSteps(1);
#endif

		// make problem and solver objects; measure time
		shared_ptr<TaylorGreenVortex2D> tgVortex = make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel);
		shared_ptr<Benchmark<2> > taylorGreen = tgVortex;
		timestart = clock();
		BenchmarkCFDSolver<2> solver(configuration, taylorGreen);
		time1 = clock() - timestart;

		try {
			solver.run();
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;
			cout << " OK ... Init: " << time1 << " sec; Run: " << time2
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
			cout << " Error" << endl;
		}

	} /* for refinement level */

	cout << "Convergence analysis (basic) terminated." << endl;

	return 0;
}
