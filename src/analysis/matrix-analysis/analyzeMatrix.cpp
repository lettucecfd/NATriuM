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

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

int main() {

	// ----------------------------------------------------------------------------------------------------

	// set Reynolds and Mach number
	const double Re = 2000;
	const double Ma = 0.1;
	const double PI = 4 * atan(1);

	// set order of FE and timeStepSize
	const size_t refinementLevel = 3;

	// ----------------------------------------------------------------------------------------------------

	//PERIODIC BOUNDARIES
	cout
			<< "Only periodic Boundaries: Calculating spectrum + pseudospectrum of streaming matrix..."
			<< endl;
	// set Problem so that the right Re and Ma are achieved
	double U = 1;
	double dqScaling = sqrt(3) * 1 / Ma;
	double viscosity = 2 * PI * U / Re; // (because L = 1)

	for (size_t orderOfFiniteElement = 2; orderOfFiniteElement <= 3;
			orderOfFiniteElement++) {

		// configure solver
		shared_ptr<SolverConfiguration> configuration = make_shared<
				SolverConfiguration>();
		std::stringstream dirname;
		dirname << getenv("NATRIUM_HOME") << "/matrix-analysis/taylorgreen"
				<< orderOfFiniteElement << "_" << refinementLevel;
		configuration->setOutputDirectory(dirname.str());
		configuration->setRestartAtLastCheckpoint(false);
		configuration->setUserInteraction(false);
		configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
		configuration->setStencilScaling(dqScaling);
		configuration->setCommandLineVerbosity(0);

		// create solver, with CFL = 1
		shared_ptr<TaylorGreenVortex2D> tgv = make_shared<TaylorGreenVortex2D>(
				viscosity, refinementLevel);
		shared_ptr<Benchmark<2> > tgBenchmark = tgv;
		configuration->setTimeStepSize(
				1.0 / dqScaling
						* CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
								*tgv->getTriangulation(), orderOfFiniteElement));
		cout << "dt = " << configuration->getTimeStepSize() << endl;
		shared_ptr<CFDSolver<2> > solver = make_shared<BenchmarkCFDSolver<2> >(
				configuration, tgBenchmark);

		// analyze eigenvalues
		matrixAnalysis<2> analyzer(solver);
		analyzer.writeSpectrum();
		analyzer.writePseudospectrum();

	}
	cout << "done." << endl;

	// ----------------------------------------------------------------------------------------------------

	//WALL BOUNDARIES
	cout
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
		shared_ptr<SolverConfiguration> configuration = make_shared<
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
		shared_ptr<CouetteFlow2D> couetteFlow = make_shared<CouetteFlow2D>(
				viscosity, U, refinementLevel, 1.0, startTime);
		shared_ptr<Benchmark<2> > couetteProblem = couetteFlow;
		configuration->setTimeStepSize(
				1.0 / dqScaling
						* CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
								*couetteFlow->getTriangulation(),
								orderOfFiniteElement));
		cout << "dt = " << configuration->getTimeStepSize() << endl;
		shared_ptr<CFDSolver<2> > solver = make_shared<BenchmarkCFDSolver<2> >(
				configuration, couetteProblem);

		// analyze eigenvalues
		matrixAnalysis<2> analyzer(solver);
		analyzer.writeSpectrum();
		analyzer.writePseudospectrum();

	}
	cout << "done." << endl;

	// ----------------------------------------------------------------------------------------------------

	return 0;
}

