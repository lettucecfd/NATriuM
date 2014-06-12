/**
 * @file SolverConfiguration_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "solver/SolverConfiguration.h"

#include <fstream>

#include "boost/test/unit_test.hpp"
#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"

#include "utilities/BasicNames.h"

namespace natrium {

#define TEST_USER_INTERACTION_CIN

BOOST_AUTO_TEST_SUITE(CFDSolverConfiguration_test)

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_Construction_test) {
	cout << "CFDSolverConfiguration_Construction_test..." << endl;
	SolverConfiguration config;
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_CreateParameterFiles_test) {
	cout << "CFDSolverConfiguration_CreateXMLFile_test..." << endl;
	SolverConfiguration config;

	std::stringstream name1;
	name1 << getenv("NATRIUM_HOME") << "/NATriuM_parameters.xml";
	std::ofstream paraOutFile1(name1.str());
	config.print_parameters(paraOutFile1, dealii::ParameterHandler::XML);

	std::stringstream name2;
	name2 << getenv("NATRIUM_HOME") << "/NATriuM_parameters.tex";
	std::ofstream paraOutFile2(name2.str());
	config.print_parameters(paraOutFile2, dealii::ParameterHandler::LaTeX);

	std::stringstream name3;
	name3 << getenv("NATRIUM_HOME") << "/NATriuM_parameters.txt";
	std::ofstream paraOutFile3(name2.str());
	config.print_parameters(paraOutFile3,
			dealii::ParameterHandler::Description);

	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_CheckSet_test) {
	cout << "CFDSolverConfiguration_CheckSet_test..." << endl;
	SolverConfiguration config;

	// Set parameter
	config.setInitializationScheme(ITERATIVE);
	BOOST_CHECK(config.getInitializationScheme() == ITERATIVE);

	cout << "done" << endl;
} /*CFDSolverConfiguration_CheckSet_test*/

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_ReadFromFile_test) {
	cout << "CFDSolverConfiguration_ReadFromFile_test..." << endl;
	SolverConfiguration config;

	try {
		config.readFromXMLFile("../src/preprocessing/NATriuM_parameters.xml");
	} catch (std::exception& e) {
		cout << e.what()
				<< "  Error! ARE YOU SURE THAT src/preprocessing/NATriuM_parameters.xml is up to date. Try to replace it by results/NATriuM_parameters.xml";
	}
	BOOST_CHECK_THROW(
			config.readFromXMLFile("../src/test/solver/invalid_parameters.xml"),
			ConfigurationException);

	cout << "done" << endl;
} /*CFDSolverConfiguration_ReadFromFile_test*/

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_PrepareOutputDirectory_test) {
	cout << "CFDSolverConfiguration_PrepareOutputDirectory_test..." << endl;
	SolverConfiguration config;

	const boost::filesystem::path natriumTmpDir("/tmp/natrium");
	const boost::filesystem::path outputDir(natriumTmpDir / "test_outputDir");

	// Sanity test
	config.setRestartAtLastCheckpoint(true);
	config.setUserInteraction(false);
	if (boost::filesystem::is_directory(natriumTmpDir)) {
		boost::filesystem::remove_all(natriumTmpDir);
	}
	boost::filesystem::create_directory(natriumTmpDir);
	config.setOutputDirectory(outputDir.string());
	config.prepareOutputDirectory();
	BOOST_CHECK(boost::filesystem::exists(outputDir));
	BOOST_CHECK_NO_THROW(config.prepareOutputDirectory());
	std::ofstream outstream((outputDir / "test.txt").string().c_str());
	outstream << "";
	BOOST_CHECK_NO_THROW(config.prepareOutputDirectory());

	// Failure test
	boost::filesystem::remove_all(natriumTmpDir);
	// Parent path not existent
	BOOST_CHECK_THROW(config.prepareOutputDirectory(), ConfigurationException);
	boost::filesystem::create_directory(natriumTmpDir);
	// No writing permissions
	boost::filesystem::perms onlyReading = boost::filesystem::others_read
			| boost::filesystem::owner_read;
	boost::filesystem::perms allPermissions = boost::filesystem::all_all;
	boost::filesystem::permissions(natriumTmpDir, onlyReading);
	// ... on parent
	BOOST_CHECK_THROW(config.prepareOutputDirectory(), ConfigurationException);
	boost::filesystem::permissions(natriumTmpDir, allPermissions);
	boost::filesystem::create_directory(outputDir);
	boost::filesystem::permissions(outputDir, onlyReading);
	// ... on output dir
	BOOST_CHECK_THROW(config.prepareOutputDirectory(), ConfigurationException);
	boost::filesystem::permissions(outputDir, allPermissions);
	boost::filesystem::create_directory(outputDir);
	std::ofstream outstream2((outputDir / "test2.txt").string().c_str());
	outstream2 << "";
	BOOST_CHECK(boost::filesystem::exists(outputDir / "test2.txt"));
	boost::filesystem::permissions(outputDir / "test2.txt", onlyReading);
	std::ofstream a((outputDir / "test2.txt").c_str());
	// .. on file in output dir
	BOOST_CHECK_THROW(config.prepareOutputDirectory(), ConfigurationException);
	boost::filesystem::permissions(outputDir, allPermissions);
	a.close();
	boost::filesystem::remove((outputDir / "test2.txt").c_str());
	// check cin
#ifdef TEST_USER_INTERACTION_CIN
	config.setRestartAtLastCheckpoint(false);
	config.setUserInteraction(true);
	config.prepareOutputDirectory();
#endif
	//clean up
	boost::filesystem::remove_all(natriumTmpDir);

	cout << "done" << endl;
} /*CFDSolverConfiguration_PrepareOutputDirectory_test*/

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_OutputFlags_test) {
	cout << "CFDSolverConfiguration_OutputFlags_test..." << endl;
	SolverConfiguration cfg;

	// Check implications
	/*cfg.setOutputFlags(out_CommandLineBasic);
	 BOOST_CHECK((out_CommandLineBasic & cfg.getOutputFlags()) != 0);
	 BOOST_CHECK((out_CommandLineError & cfg.getOutputFlags()) != 0);
	 BOOST_CHECK(not ((out_CommandLineFull & cfg.getOutputFlags()) != 0));

	 cfg.setOutputFlags(out_CommandLineFull);
	 BOOST_CHECK((out_CommandLineBasic & cfg.getOutputFlags()) != 0);
	 BOOST_CHECK((out_CommandLineError & cfg.getOutputFlags()) != 0);
	 BOOST_CHECK((out_CommandLineFull & cfg.getOutputFlags()) != 0);
	 */
	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
