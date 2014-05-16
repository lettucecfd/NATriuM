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

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_Construction_test){
	cout << "CFDSolverConfiguration_Construction_test..." << endl;
	SolverConfiguration config;
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_CreateXMLFile_test){
	cout << "CFDSolverConfiguration_CreateXMLFile_test..." << endl;
	SolverConfiguration config;

	std::ofstream paraOutFile("../results/NATriuM_parameters.xml");
	config.print_parameters(paraOutFile, dealii::ParameterHandler::XML);
	cout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(CFDSolverConfiguration_OutputFlags_test){
	cout << "CFDSolverConfiguration_OutputFlags_test..." << endl;
	SolverConfiguration cfg;

	// Check implications
	cfg.setOutputFlags(out_CommandLineBasic);
	BOOST_CHECK((out_CommandLineBasic & cfg.getOutputFlags()) != 0);
	BOOST_CHECK((out_CommandLineError & cfg.getOutputFlags()) != 0);
	BOOST_CHECK(not ((out_CommandLineFull & cfg.getOutputFlags()) != 0));

	cfg.setOutputFlags(out_CommandLineFull);
	BOOST_CHECK((out_CommandLineBasic & cfg.getOutputFlags()) != 0);
	BOOST_CHECK((out_CommandLineError & cfg.getOutputFlags()) != 0);
	BOOST_CHECK((out_CommandLineFull & cfg.getOutputFlags()) != 0);

	cout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
