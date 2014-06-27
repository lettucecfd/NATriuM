#include "solver/Stats.h"

#include "boost/test/unit_test.hpp"
#include "boost/filesystem.hpp"

#include "utilities/BasicNames.h"

#include "PeriodicTestFlow2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(Stats_test)

BOOST_AUTO_TEST_CASE(SolverStats_ConstructionAndFunctionality_test) {
	cout << "SolverStats_ConstructionAndFunctionality_test..." << endl;

	// make solver object
	const boost::filesystem::path natriumTmpDir("/tmp/natrium_test");
	const boost::filesystem::path natriumTmpFile = natriumTmpDir / "testtable.txt";
	shared_ptr<SolverConfiguration> testConfiguration = make_shared<SolverConfiguration>();
	testConfiguration->setOutputDirectory(natriumTmpDir.c_str());
	testConfiguration->setCommandLineVerbosity(SILENT);
	testConfiguration->setUserInteraction(false);
	testConfiguration->setNumberOfTimeSteps(10);
	size_t refinementLevel = 3;
	double viscosity = 0.9;
	shared_ptr<ProblemDescription<2> > testFlow = make_shared<SteadyPeriodicTestFlow2D>(viscosity, refinementLevel);
	CFDSolver<2> solver(testConfiguration, testFlow);

	// make solver stats object
	if (boost::filesystem::is_directory(natriumTmpDir)) {
		boost::filesystem::remove_all(natriumTmpDir);
	}
	boost::filesystem::create_directory(natriumTmpDir);
	SolverStats<2> stats(solver, natriumTmpFile.c_str());
	BOOST_CHECK(boost::filesystem::exists(natriumTmpFile));

	// check uptodate
	BOOST_CHECK(not stats.isUpToDate());
	stats.update();
	BOOST_CHECK(stats.isUpToDate());

	// count lines
	std::ifstream file(natriumTmpFile.c_str());
	size_t linecount = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	BOOST_CHECK(linecount == 1);

	// check printNewLine
	solver.run();
	stats.printNewLine();
	BOOST_CHECK_EQUAL(stats.getIterationNumber(), testConfiguration->getNumberOfTimeSteps());
	linecount = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
	BOOST_CHECK(linecount > 1);


	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
