/**
 * @file SolverConfiguration_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/solver/SolverConfiguration.h"

#include <fstream>

#include "boost/test/unit_test.hpp"
#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

//#define TEST_USER_INTERACTION_CIN

BOOST_AUTO_TEST_SUITE(CFDSolverConfiguration_test)

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_Construction_test) {
	pout << "CFDSolverConfiguration_Construction_test..." << endl;
	SolverConfiguration config;
	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_CreateParameterFiles_test) {
	pout << "CFDSolverConfiguration_CreateXMLFile_test..." << endl;
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

	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_CheckSet_test) {
	pout << "CFDSolverConfiguration_CheckSet_test..." << endl;
	SolverConfiguration config;

	// Set parameter
	pout << " ... sanity test ... " << endl;
	config.setInitializationScheme(ITERATIVE);
	BOOST_CHECK(config.getInitializationScheme() == ITERATIVE);
	config.setIterativeInitializationNumberOfIterations(10);
	BOOST_CHECK_EQUAL(config.getIterativeInitializationNumberOfIterations(),
			(size_t ) 10);
	config.setIterativeInitializationResidual(0.001);
	BOOST_CHECK(0.001 == config.getIterativeInitializationResidual());

	config.setSimulationEndTime(10.0);
	BOOST_CHECK_CLOSE(config.getSimulationEndTime(), 10.0, 1e-10);
	config.setNumberOfTimeSteps(100);
	BOOST_CHECK_EQUAL(config.getNumberOfTimeSteps(), (size_t ) 100);
	config.setConvergenceThreshold(1e-12);
	BOOST_CHECK_CLOSE(config.getConvergenceThreshold(), 1e-12, 1e-5);
	config.setBGKSteadyStateGamma(0.1);
	BOOST_CHECK_CLOSE(config.getBGKSteadyStateGamma(), 0.1, 1e-5);
	config.setCollisionScheme(BGK_STEADY_STATE);
	BOOST_CHECK_EQUAL(config.getCollisionScheme(), BGK_STEADY_STATE);
	BOOST_CHECK_EQUAL(config.getForcingScheme(), NO_FORCING);
	config.setForcingScheme(SHIFTING_VELOCITY);
	BOOST_CHECK_EQUAL(config.getForcingScheme(), SHIFTING_VELOCITY);
	config.setForcingScheme(EXACT_DIFFERENCE);
	BOOST_CHECK_EQUAL(config.getForcingScheme(), EXACT_DIFFERENCE);
	config.setForcingScheme(GUO);
	BOOST_CHECK_EQUAL(config.getForcingScheme(), GUO);
	BOOST_CHECK_EQUAL(config.isFiltering(), false);
	config.setFiltering(true);
	BOOST_CHECK_EQUAL(config.isFiltering(), true);
	BOOST_CHECK_EQUAL(config.getFilteringScheme(), EXPONENTIAL_FILTER);
	config.setFilteringScheme(NEW_FILTER);
	BOOST_CHECK_EQUAL(config.getFilteringScheme(), NEW_FILTER);
	BOOST_CHECK_EQUAL(config.getExponentialFilterAlpha(), 10.0);
	config.setExponentialFilterAlpha(2.0);
	BOOST_CHECK_EQUAL(config.getExponentialFilterAlpha(), 2.0);
	BOOST_CHECK_EQUAL(config.getExponentialFilterS(), 20.0);
	config.setExponentialFilterS(2.0);
	BOOST_CHECK_EQUAL(config.getExponentialFilterS(), 2.0);
	BOOST_CHECK_EQUAL(config.isOutputTurbulenceStatistics(), false);
	config.setOutputTurbulenceStatistics(true);
	BOOST_CHECK_EQUAL(config.isOutputTurbulenceStatistics(), true);
	BOOST_CHECK_EQUAL(config.getWallNormalDirection(), size_t(1));
	config.setWallNormalDirection(0);
	BOOST_CHECK_EQUAL(config.getWallNormalDirection(), size_t(0));
	vector<double> coord = config.getWallNormalCoordinates();
	BOOST_CHECK_EQUAL(coord.size(), size_t(3));
	BOOST_CHECK_CLOSE(coord.at(0), 0.1, 1e-4);
	BOOST_CHECK_CLOSE(coord.at(1), 0.2, 1e-4);
	BOOST_CHECK_CLOSE(coord.at(2), 0.5, 1e-4);
	coord.at(0) = 0.05;
	config.setWallNormalCoordinates(coord);
	vector<double> coord2 = config.getWallNormalCoordinates();
	BOOST_CHECK_CLOSE(coord2.at(0), 0.05, 1e-8);
	BOOST_CHECK_CLOSE(config.getCFL(), 0.4, 1e-10);
	config.setCFL(10);
	BOOST_CHECK_CLOSE(config.getCFL(), 10, 1e-10);
	BOOST_CHECK_EQUAL(config.getRestartAtIteration(), size_t(0));
	config.setRestartAtIteration(1);
	BOOST_CHECK_EQUAL(config.getRestartAtIteration(), size_t(1));
	BOOST_CHECK_EQUAL(config.getAdvectionScheme(), SEDG);
	config.setAdvectionScheme(SEMI_LAGRANGIAN);
	BOOST_CHECK_EQUAL(config.getAdvectionScheme(), SEMI_LAGRANGIAN);
	BOOST_CHECK_EQUAL(config.getFilterInterval(), size_t(1));
	config.setFilterInterval(2);
	BOOST_CHECK_EQUAL(config.getFilterInterval(), size_t(2));
	BOOST_CHECK_EQUAL(config.isFilterDegreeByComponentSums(), false);
	config.setFilterDegreeByComponentSums(true);
	BOOST_CHECK_EQUAL(config.isFilterDegreeByComponentSums(), true);
	BOOST_CHECK_EQUAL(config.isOutputGlobalTurbulenceStatistics(), false);
	config.setOutputGlobalTurbulenceStatistics(true);
	BOOST_CHECK_EQUAL(config.isOutputGlobalTurbulenceStatistics(), true);
	BOOST_CHECK_EQUAL(config.isVmultLimiter(), false);
	config.setVmultLimiter(true);
	BOOST_CHECK_EQUAL(config.isVmultLimiter(), true);
	BOOST_CHECK_EQUAL(config.getQuadrature(), QGAUSS_LOBATTO);
	config.setQuadrature(QGAUSS);
	BOOST_CHECK_EQUAL(config.getQuadrature(), QGAUSS);
	BOOST_CHECK_EQUAL(config.getSupportPoints(), GAUSS_LOBATTO_POINTS);
	config.setSupportPoints(EQUIDISTANT_POINTS);
	BOOST_CHECK_EQUAL(config.getSupportPoints(), EQUIDISTANT_POINTS);
	config.setSupportPoints(GAUSS_CHEBYSHEV_POINTS);
	BOOST_CHECK_EQUAL(config.getSupportPoints(), GAUSS_CHEBYSHEV_POINTS);
	config.setSupportPoints(GAUSS_LOBATTO_CHEBYSHEV_POINTS);
	BOOST_CHECK_EQUAL(config.getSupportPoints(), GAUSS_LOBATTO_CHEBYSHEV_POINTS);
	BOOST_CHECK_EQUAL(config.getMRTBasis(), DELLAR_D2Q9);
	config.setMRTBasis(LALLEMAND_D2Q9);
	BOOST_CHECK_EQUAL(config.getMRTBasis(), LALLEMAND_D2Q9);
	config.setMRTBasis(DHUMIERES_D3Q19);
	BOOST_CHECK_EQUAL(config.getMRTBasis(), DHUMIERES_D3Q19);
	BOOST_CHECK_EQUAL(config.getMRTRelaxationTimes(), RELAX_FULL);
	config.setMRTRelaxationTimes(DELLAR_RELAX_ONLY_N);
	BOOST_CHECK_EQUAL(config.getMRTRelaxationTimes(), DELLAR_RELAX_ONLY_N);
	config.setMRTRelaxationTimes(RELAX_DHUMIERES_PAPER);
	BOOST_CHECK_EQUAL(config.getMRTRelaxationTimes(), RELAX_DHUMIERES_PAPER);
	BOOST_CHECK_EQUAL(config.getEquilibriumScheme(), BGK_EQUILIBRIUM);
	config.setEquilibriumScheme(ENTROPIC_EQUILIBRIUM);
	BOOST_CHECK_EQUAL(config.getEquilibriumScheme(), ENTROPIC_EQUILIBRIUM);
	config.setEquilibriumScheme(INCOMPRESSIBLE_EQUILIBRIUM);
	BOOST_CHECK_EQUAL(config.getEquilibriumScheme(), INCOMPRESSIBLE_EQUILIBRIUM);
	config.setEquilibriumScheme(STEADYSTATE_EQUILIBRIUM);
	BOOST_CHECK_EQUAL(config.getEquilibriumScheme(), STEADYSTATE_EQUILIBRIUM);

	/// Failure test
	pout << " ... failure test ... " << endl;
	BOOST_CHECK_THROW(config.setSimulationEndTime(-0.1),
			ConfigurationException);
	BOOST_CHECK_THROW(config.setNumberOfTimeSteps(-0.1),
			ConfigurationException);
	BOOST_CHECK_THROW(config.setBGKSteadyStateGamma(-0.1),
			ConfigurationException);
	BOOST_CHECK_THROW(config.setCFL(-1.0), ConfigurationException);
	BOOST_CHECK_THROW(config.setRestartAtIteration(-10),
			ConfigurationException);

	pout << "done" << endl;
} /*CFDSolverConfiguration_CheckSet_test*/

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_ReadFromFile_test) {

	pout
			<< "CFDSolverConfiguration_CheckIfEnvironmentVariable_NATRIUM_DIR_IsSet_text"
			<< endl;
	std::stringstream name1;
	name1 << getenv("NATRIUM_DIR");
	BOOST_CHECK(name1.str().size() != 0);
	pout << "done." << endl;

	pout << "CFDSolverConfiguration_ReadFromFile_test..." << endl;
	SolverConfiguration config;
	name1 << "/src/preprocessing/NATriuM_parameters.xml";
	std::stringstream name2;
	name2 << getenv("NATRIUM_DIR");
	name2 << "/src/test/solver/invalid_parameters.xml";

	try {
		config.readFromXMLFile(name1.str());
	} catch (std::exception& e) {
		pout << e.what()
				<< "  Error! ARE YOU SURE THAT src/preprocessing/NATriuM_parameters.xml is up to date? "
						"Try to replace it by results/NATriuM_parameters.xml"
				<< endl;
	}
	BOOST_CHECK_THROW(config.readFromXMLFile(name2.str()),
			ConfigurationException);

	pout << "done" << endl;
} /*CFDSolverConfiguration_ReadFromFile_test*/

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_PrepareOutputDirectory_test) {
	pout << "CFDSolverConfiguration_PrepareOutputDirectory_test..." << endl;
	SolverConfiguration config;

	const boost::filesystem::path natriumTmpDir("/tmp/natrium");
	const boost::filesystem::path outputDir(natriumTmpDir / "test_outputDir");

	// Sanity test
	config.setRestartAtIteration(1);
	config.setUserInteraction(false);
	if (is_MPI_rank_0()){
		if (boost::filesystem::is_directory(natriumTmpDir)) {
			boost::filesystem::remove_all(natriumTmpDir);
		}
		boost::filesystem::create_directory(natriumTmpDir);
	}

	// wait for directory to be created
	sleep(1);
	MPI_sync();
	config.setOutputDirectory(outputDir.string());
	config.prepareOutputDirectory();
	MPI_sync();
	BOOST_CHECK(boost::filesystem::exists(outputDir));
	BOOST_CHECK_NO_THROW(config.prepareOutputDirectory());
	if (is_MPI_rank_0()) {
		std::ofstream outstream((outputDir / "test.txt").string().c_str());
		outstream << "";
		outstream.close();
	}
	config.prepareOutputDirectory();
	MPI_sync();
	BOOST_CHECK_NO_THROW(config.prepareOutputDirectory());

	// Failure test (failures are only detected in rank 0)
	if (is_MPI_rank_0()) {
		boost::filesystem::remove_all(natriumTmpDir);
		// Parent path not existent
		BOOST_CHECK_THROW(config.prepareOutputDirectory(),
				ConfigurationException);
		boost::filesystem::create_directory(natriumTmpDir);
		// No writing permissions
		boost::filesystem::perms onlyReading = boost::filesystem::others_read
				| boost::filesystem::owner_read;
		boost::filesystem::perms allPermissions = boost::filesystem::all_all;
		boost::filesystem::permissions(natriumTmpDir, onlyReading);
		// ... on parent
		BOOST_CHECK_THROW(config.prepareOutputDirectory(),
				ConfigurationException);
		boost::filesystem::permissions(natriumTmpDir, allPermissions);
		boost::filesystem::create_directory(outputDir);
		boost::filesystem::permissions(outputDir, onlyReading);

		// ... on output dir
		BOOST_CHECK_THROW(config.prepareOutputDirectory(),
				ConfigurationException);
		boost::filesystem::permissions(outputDir, allPermissions);
		boost::filesystem::create_directory(outputDir);
		std::ofstream outstream2((outputDir / "test2.txt").string().c_str());
		outstream2 << "";
		outstream2.close();
		BOOST_CHECK(boost::filesystem::exists(outputDir / "test2.txt"));
		boost::filesystem::permissions(outputDir / "test2.txt", onlyReading);
		// .. on file in output dir
		std::ofstream a((outputDir / "test2.txt").c_str());
		BOOST_CHECK_THROW(config.prepareOutputDirectory(),
				ConfigurationException);
		boost::filesystem::permissions(outputDir, allPermissions);
		a.close();
		boost::filesystem::remove((outputDir / "test2.txt").c_str());
		// check cin
#ifdef TEST_USER_INTERACTION_CIN
		config.setRestartAtIteration(0);
		config.setUserInteraction(true);
		config.prepareOutputDirectory();
#endif
		//clean up
		boost::filesystem::remove_all(natriumTmpDir);
	}
	pout << "done" << endl;
} /*CFDSolverConfiguration_PrepareOutputDirectory_test*/

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_OutputFlags_test) {
	pout << "CFDSolverConfiguration_OutputFlags_test..." << endl;
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
	pout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

