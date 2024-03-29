/*
 * ErrorStats_test.cpp
 *
 *  Created on: 27.06.2014
 *      Author: kraemer
 */

#include "natrium/solver/ErrorStats.h"
#include <fstream>

#include "boost/test/unit_test.hpp"
#include "boost/filesystem.hpp"

#include "natrium/solver/SolverConfiguration.h"
#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/problemdescription/Benchmark.h"
#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/TaylorGreenVortex2D.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(ErrorStats_test)

BOOST_AUTO_TEST_CASE(ErrorStats_ConstructionAndFunctionality_test) {
	pout << "ErrorStats_ConstructionAndFunctionality_test..." << endl;

	// make solver object
	const boost::filesystem::path natriumTmpDir("/tmp/natrium_test");
	const boost::filesystem::path natriumTmpFile = natriumTmpDir
			/ "test_errorstable.txt";
	boost::shared_ptr<SolverConfiguration> testConfiguration = boost::make_shared<
			SolverConfiguration>();
	testConfiguration->setOutputDirectory(natriumTmpDir.c_str());
	testConfiguration->setCommandLineVerbosity(SILENT);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setNumberOfTimeSteps(10);
	size_t refinementLevel = 3;
	double viscosity = 0.9;
	boost::shared_ptr<Benchmark<2> > testFlow = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel);
	BenchmarkCFDSolver<2> solver(testConfiguration, testFlow);

	// make error stats object
	if (is_MPI_rank_0()) {
		if (boost::filesystem::is_directory(natriumTmpDir)) {
			boost::filesystem::remove_all(natriumTmpDir);
		}
		boost::filesystem::create_directory(natriumTmpDir);
		boost::filesystem::create_directory(natriumTmpDir / "checkpoint");
	}
	ErrorStats<2> stats(&solver, natriumTmpFile.c_str());
	BOOST_CHECK(boost::filesystem::exists(natriumTmpFile));

	// check uptodate
	BOOST_CHECK(not stats.isUpToDate());
	stats.update();
	BOOST_CHECK(stats.isUpToDate());

	// count lines
	std::ifstream file(natriumTmpFile.c_str());
	MPI_sync();
	size_t linecount = std::count(std::istreambuf_iterator<char>(file),
			std::istreambuf_iterator<char>(), '\n');
	BOOST_CHECK(linecount == 1);

	// check printNewLine
	stats.printNewLine();
	solver.run();
	stats.printNewLine();
	MPI_sync();
	BOOST_CHECK_EQUAL(stats.getIterationNumber(),
			testConfiguration->getNumberOfTimeSteps());
	file.seekg(0); //set read position to 0
	linecount = std::count(std::istreambuf_iterator<char>(file),
			std::istreambuf_iterator<char>(), '\n');
	BOOST_CHECK_GT(linecount, 1);

	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(ErrorStats_FunctionsInSolverContext_test) {
	pout << "ErrorStats_FunctionsInSolverContext_test..." << endl;
	// make solver object
	const boost::filesystem::path natriumTmpDir("/tmp/natrium_test");
	boost::shared_ptr<SolverConfiguration> testConfiguration = boost::make_shared<
			SolverConfiguration>();
	testConfiguration->setOutputDirectory(natriumTmpDir.c_str());
	testConfiguration->setCommandLineVerbosity(SILENT);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setNumberOfTimeSteps(10);
	testConfiguration->setOutputTableInterval(1);
	size_t refinementLevel = 3;
	double viscosity = 0.9;
	boost::shared_ptr<Benchmark<2> > testFlow = boost::make_shared<TaylorGreenVortex2D>(
			viscosity, refinementLevel);
	BenchmarkCFDSolver<2> solver(testConfiguration, testFlow);

	solver.run();
	solver.getSolverStats()->update();
	BOOST_CHECK_EQUAL(solver.getSolverStats()->getIterationNumber(),
			testConfiguration->getNumberOfTimeSteps());

	// count lines
	std::ifstream file(solver.getErrorStats()->getFilename());
	MPI_sync();
	size_t linecount = std::count(std::istreambuf_iterator<char>(file),
			std::istreambuf_iterator<char>(), '\n');
	BOOST_CHECK(linecount > 8);

	pout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

