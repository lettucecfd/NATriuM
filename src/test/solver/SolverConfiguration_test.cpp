/**
 * @file SolverConfiguration_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "solver/SolverConfiguration.h"

#include <fstream>

#include "boost/test/unit_test.hpp"

#include "utilities/BasicNames.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(CFDSolverConfiguration_test)

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_Construction_test) {
	cout << "CFDSolverConfiguration_Construction_test..." << endl;
	SolverConfiguration config;
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_CreateParameterFiles_test) {
	cout << "CFDSolverConfiguration_CreateXMLFile_test..." << endl;
	SolverConfiguration config;

	std::ofstream paraOutFile1("../results/NATriuM_parameters.xml");
	config.print_parameters(paraOutFile1, dealii::ParameterHandler::XML);

	std::ofstream paraOutFile2("../results/NATriuM_parameters.tex");
	config.print_parameters(paraOutFile2, dealii::ParameterHandler::LaTeX);

	std::ofstream paraOutFile3("../results/NATriuM_parameters.txt");
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
		cout << e.what() << "  Error! ARE YOU SURE THAT src/preprocessing/NATriuM_parameters.xml is up to date. Try to replace it by results/NATriuM_parameters.xml";
	}
	BOOST_CHECK_THROW(
			config.readFromXMLFile("../src/test/solver/invalid_parameters.xml"),
			ConfigurationException);

	cout << "done" << endl;
} /*CFDSolverConfiguration_ReadFromFile_test*/

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
