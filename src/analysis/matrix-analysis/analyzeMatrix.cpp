/**
 * @file analyzeMatrix.cpp
 * @short Executable for the analysis of spectra/pseudospectra of streaming matrices
 * @date 01.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "matrixAnalysis.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/problemdescription/Benchmark.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"
#include "natrium/benchmarks/CouetteFlow2D.h"

#include "natrium/stencils/D2Q9.h"

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

int main() {

	MPIGuard::getInstance();

	// ----------------------------------------------------------------------------------------------------

	const bool semi_lagrange = 1;

	// set Reynolds and Mach number
	const double Re = 2000;
	const double Ma = 0.1;
	const double PI = 4 * atan(1);

	// set order of FE and timeStepSize
	const size_t refinementLevel = 2;

	// ----------------------------------------------------------------------------------------------------

	//PERIODIC BOUNDARIES
	pout
			<< "Only periodic Boundaries: Calculating spectrum + pseudospectrum of streaming matrix..."
			<< endl;
	// set Problem so that the right Re and Ma are achieved
	double U = 1;
	double dqScaling = sqrt(3) * 1 / Ma;
	double viscosity = 2 * PI * U / Re; // (because L = 1)

	for (size_t orderOfFiniteElement = 1; orderOfFiniteElement <= 9;
			orderOfFiniteElement++) {

		// configure solver
		boost::shared_ptr<SolverConfiguration> configuration =
				boost::make_shared<SolverConfiguration>();
		std::stringstream dirname;
		dirname << getenv("NATRIUM_HOME") << "/matrix-analysis/taylorgreen"
				<< orderOfFiniteElement << "_" << refinementLevel;
		configuration->setOutputDirectory(dirname.str());
		configuration->setUserInteraction(false);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(dqScaling);
		configuration->setCommandLineVerbosity(0);
		configuration->setCFL(10);
		if (semi_lagrange) {
			configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
		}

		// create solver, with CFL = 1
		boost::shared_ptr<TaylorGreenVortex2D> tgv = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, refinementLevel);
		boost::shared_ptr<Benchmark<2> > tgBenchmark = tgv;
		boost::shared_ptr<CFDSolver<2> > solver = boost::make_shared<
				BenchmarkCFDSolver<2> >(configuration, tgBenchmark);
		pout << "dt = " << solver->getTimeStepSize() << endl;

		// analyze eigenvalues
		matrixAnalysis<2> analyzer(solver);
		vector<std::complex<double> > eigenvalues;
		if (semi_lagrange) {
			pout << "abs max: "
					<< analyzer.computeSpectrum(
									solver->getAdvectionOperator()->getSystemMatrix(),
									eigenvalues, 0.0) << endl;
			analyzer.writeSpectrum(false);
			analyzer.writePseudospectrum(false, 10, 0.01);
		} else {
			pout << "abs max: "
					<< solver->getTimeStepSize()
							* analyzer.computeSpectrum(
									solver->getAdvectionOperator()->getSystemMatrix(),
									eigenvalues, 0.0) << endl;
			analyzer.writeSpectrum(true);
			analyzer.writePseudospectrum(true);
		}

	}
	pout << "done." << endl;

	// ----------------------------------------------------------------------------------------------------
	/*
	 //WALL BOUNDARIES
	 pout
	 << "With wall Boundaries: Calculating spectrum + pseudospectrum of streaming matrix..."
	 << endl;

	 // set Problem so that the right Re and Ma are achieved
	 // set consistent U:     U/L (couette) = U/L (taylor-green) = 1/2PI
	 U = 1 / (2 * PI);	// / sqrt(3) * Ma;
	 dqScaling = sqrt(3) * U / Ma;
	 viscosity = U / Re; // (because L = 1)
	 // in order to start from a continuous solution, do not start at t=0
	 const double startTime = 0.0;
	 // set small time step size

	 for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 3;
	 orderOfFiniteElement++) {

	 // configure solver
	 boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
	 SolverConfiguration>();
	 std::stringstream dirname;
	 dirname << getenv("NATRIUM_HOME") << "/matrix-analysis/couette"
	 << orderOfFiniteElement << "_" << refinementLevel;
	 configuration->setOutputDirectory(dirname.str());
	 configuration->setRestartAtLastCheckpoint(false);
	 configuration->setUserInteraction(false);
	 configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	 configuration->setStencilScaling(dqScaling);
	 configuration->setCommandLineVerbosity(0);

	 // create solver, with CFL = 1
	 boost::shared_ptr<CouetteFlow2D> couetteFlow = boost::make_shared<CouetteFlow2D>(
	 viscosity, U, refinementLevel, 1.0, startTime);
	 boost::shared_ptr<Benchmark<2> > couetteProblem = couetteFlow;
	 configuration->setTimeStepSize(
	 1.0 / dqScaling
	 * CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
	 *couetteFlow->getMesh(),
	 orderOfFiniteElement));
	 pout << "dt = " << configuration->getTimeStepSize() << endl;
	 boost::shared_ptr<CFDSolver<2> > solver = boost::make_shared<BenchmarkCFDSolver<2> >(
	 configuration, couetteProblem);

	 // analyze eigenvalues
	 matrixAnalysis<2> analyzer(solver);
	 analyzer.writeSpectrum();
	 analyzer.writePseudospectrum();

	 }
	 pout << "done." << endl;
	 */
	// ----------------------------------------------------------------------------------------------------
	return 0;
}

