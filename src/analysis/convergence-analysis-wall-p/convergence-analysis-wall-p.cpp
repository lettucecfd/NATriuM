/**
 * @file convergence-analysis-p.cpp
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

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/CouetteFlow2D.h"

using namespace natrium;

// if this define statement is enabled: only the initialization time is regarded
//#define MEASURE_ONLY_INIT_TIME

// Main function
int main() {

	MPIGuard::getInstance();

	pout
			<< "Starting NATriuM convergence analysis with moving wall boundaries (various p) ..."
			<< endl;

	/////////////////////////////////////////////////
	// set parameters, set up configuration object
	//////////////////////////////////////////////////

	// specify Reynolds number
	const double Re = 2000;
	// specify Mach number
	const double Ma = 0.05;

	// Problem description
	// -------------------
	double scaling = 1; //sqrt(3) * U / Ma;
	// velocity of top plate
	const double U = scaling * Ma / sqrt(3); //0.01;
	// length of quadratic domain
	const double L = 1.0;
	// Viscosity
	const double viscosity = U * L / Re;
	// starting time
	//const double t0 = 30.0;
	const double t0 = 0.0;

	size_t refinementLevel = 3;
	// prepare time table file
	// the output is written to the standard output directory (e.g. NATriuM/results or similar)
	std::stringstream filename;
	filename << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-wall-p/table_runtime.txt";
	std::ofstream timeFile(filename.str().c_str());
	timeFile
			<< "# order of FE   dt        init time (sec)             loop time (sec)         time for one iteration (sec)"
			<< endl;

	// prepare error table file
	std::stringstream filename2;
	filename2 << getenv("NATRIUM_HOME")
			<< "/convergence-analysis-wall-p/table_order.txt";
	std::ofstream orderFile(filename2.str().c_str());
	orderFile << "# visc = " << viscosity << "; Ma = " << Ma << endl;
	orderFile
			<< "#  orderOfFe  i      t         max |u_analytic|  max |error_u|  max |error_rho|   ||error_u||_2   ||error_rho||_2"
			<< endl;

	for (size_t orderOfFiniteElement = 4; orderOfFiniteElement <= 12;
			orderOfFiniteElement += 2) {
		timeFile << endl << "# order FE = " << orderOfFiniteElement << endl;
		orderFile << endl << "# order FE = " << orderOfFiniteElement << endl;
		pout << endl << "order FE = " << orderOfFiniteElement << endl;

		double dt = 0.0005;
		//for (double dt = 0.1; dt >= 0.00001; dt /= 2.) {
		timeFile << "# dt = " << dt << endl;
		orderFile << "# dt = " << dt << endl;
		pout << "dt = " << dt << endl;

		// time measurement variables
		double time1, time2, timestart;

		// setup configuration
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME") << "/convergence-analysis-wall-p/"
				<< orderOfFiniteElement << "_" << refinementLevel << "_" << dt;
		boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
				SolverConfiguration>();
		//configuration->setSwitchOutputOff(true);
		configuration->setOutputDirectory(dirName.str());
		//configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setOutputTableInterval(1000);
		//configuration->setOutputCheckpointInterval(1000);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(scaling);
		configuration->setCommandLineVerbosity(0);
		// DO NOT USE!!!!
		configuration->setCFL(0.4);
		configuration->setNumberOfTimeSteps(40.0 / dt);

#ifdef MEASURE_ONLY_INIT_TIME
		configuration->setNumberOfTimeSteps(1);
#endif

		// make problem and solver objects; measure time
		boost::shared_ptr<CouetteFlow2D> couette2D = boost::make_shared<CouetteFlow2D>(
				viscosity, U, refinementLevel, L, t0);
		boost::shared_ptr<Benchmark<2> > benchmark = couette2D;
		timestart = clock();
		BenchmarkCFDSolver<2> solver(configuration, benchmark);
		time1 = clock() - timestart;

		try {
			solver.run();
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;
			pout << " OK ... Init: " << time1 << " sec; Run: " << time2
					<< " sec." << endl;
			// put out runtime
			timeFile << orderOfFiniteElement << "         " << dt << "      "
					<< time1 << "     " << time2 << "        "
					<< time2 / configuration->getNumberOfTimeSteps() << endl;
			// put out final errors
			solver.getErrorStats()->update();
			orderFile << orderOfFiniteElement << " " << solver.getIteration()
					<< " " << solver.getTime() << " "
					<< solver.getErrorStats()->getMaxUAnalytic() << " "
					<< solver.getErrorStats()->getMaxVelocityError() << " "
					<< solver.getErrorStats()->getMaxDensityError() << " "
					<< solver.getErrorStats()->getL2VelocityError() << " "
					<< solver.getErrorStats()->getL2DensityError() << endl;
		} catch (std::exception& e) {
			pout << " Error" << endl;
		}

		//} /* for time step*/
		timeFile << endl;
		orderFile << endl;
	} /* for refinement level */

	pout << "Convergence analysis with moving walls terminated." << endl;

	return 0;
}
