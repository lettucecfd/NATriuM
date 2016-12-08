/**
 * @file CommandLineParser_test.cpp
 * @short
 * @date 06.12.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/CommandLineParser.h"

#include "boost/test/unit_test.hpp"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(CommandLineParser_test)

BOOST_AUTO_TEST_CASE(CommandLineParser_Construction_test) {
	pout << "CommandLineParser_Construction_test..." << endl;

	char arg1[] = "program_name";
	char* argv[] = { arg1, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	pout << "done." << endl;
} /* CommandLineParser_Construction_test */

BOOST_AUTO_TEST_CASE(CommandLineParser_DefineParameters_test) {
	pout << "CommandLineParser_DefineParameters_test..." << endl;

	char arg1[] = "program_name";
	char arg2[] = "--top";
	char arg3[] = "2";
	char arg4[] = "--omma";
	char arg5[] = "omi";
	char arg6[] = "--mogli";
	char arg7[] = "1.0";
	char* argv[] = { arg1, arg2, arg3, arg4, arg5, arg6, arg7, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.setArgument<int>("top", "test int", 2);
	c.setArgument<string>("omma", "test string", "");
	c.setArgument<double>("mogli", "test double", 1.0);
	c.setArgument<int>("def", "test default", 12);

	pout << "done." << endl;
} /* CommandLineParser_DefineParameters_test */

BOOST_AUTO_TEST_CASE(CommandLineParser_Help_test) {
	pout << "CommandLineParser_Help_test..." << endl;

	char arg1[] = "program_name";
	char arg6[] = "--help";
	char* argv[] = { arg1, arg6, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.setArgument<double>("double1", "test double", 1.0);
	c.setPositionalArgument<string>("string1", "test string");
	c.addDocumentationString("test_program", "a test program");

	BOOST_CHECK_THROW(c.importOptions(), HelpMessageStop);

	pout << "done." << endl;
} /* CommandLineParser_ParseDouble_test */

BOOST_AUTO_TEST_CASE(CommandLineParser_ParseDouble_test) {
	pout << "CommandLineParser_ParseDouble_test..." << endl;
	{
		char arg1[] = "program_name";
		char arg6[] = "--mogli";
		char arg7[] = "1.0";
		char* argv[] = { arg1, arg6, arg7, NULL };
		int argc = sizeof(argv) / sizeof(char*) - 1;

		CommandLineParser c(argc, argv);

		// define arguments
		c.setArgument<double>("mogli", "test double", 1.0);
		c.setArgument<double>("def", "test default", 12);
		c.importOptions();

		// read specified arguments
		BOOST_CHECK_EQUAL(1.0, c.getArgument<double>("mogli"));
		// read default arguments
		BOOST_CHECK_EQUAL(12, c.getArgument<double>("def"));

		pout << "done." << endl;
		pout << "CommandLineParser_Destruction_test..." << endl;
	}
	pout << "done." << endl;
} /* CommandLineParser_ParseDouble_test */

BOOST_AUTO_TEST_CASE(CommandLineParser_ParseInt_test) {
	pout << "CommandLineParser_ParseInt_test..." << endl;

	char arg1[] = "program_name";
	char arg2[] = "--top";
	char arg3[] = "2";
	char* argv[] = { arg1, arg2, arg3, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.setArgument<int>("top", "test int", 13);
	//c.setStringArgument("omma", "test string", "");
	//c.setDoubleArgument("mogli", "test double", 1.0);
	c.setArgument<int>("def", "test default", 12);
	c.importOptions();

	BOOST_CHECK_EQUAL(2, c.getArgument<int>("top"));
	//BOOST_CHECK_EQUAL("omi", c.getStringArgument("omma"));
	//BOOST_CHECK_EQUAL(1.0, c.getDoubleArgument("mogli"));
	BOOST_CHECK_EQUAL(12, c.getArgument<int>("def"));

	pout << "done." << endl;
} /* CommandLineParser_ParseInt_test */

BOOST_AUTO_TEST_CASE(CommandLineParser_ParseString_test) {
	pout << "CommandLineParser_ParseString_test..." << endl;

	char arg1[] = "program_name";
	char arg2[] = "--top";
	char arg3[] = "2";
	char* argv[] = { arg1, arg2, arg3, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.setArgument<string>("top", "test int", "default_string");
	//c.setStringArgument("omma", "test string", "");
	//c.setDoubleArgument("mogli", "test double", 1.0);
	c.setArgument<string>("def", "test default", "default_string2");
	c.importOptions();

	BOOST_CHECK_EQUAL("2", c.getArgument<string>("top"));
	//BOOST_CHECK_EQUAL("omi", c.getStringArgument("omma"));
	//BOOST_CHECK_EQUAL(1.0, c.getDoubleArgument("mogli"));
	BOOST_CHECK_EQUAL("default_string2", c.getArgument<string>("def"));

	pout << "done." << endl;
} /* CommandLineParser_ParseString_test */


BOOST_AUTO_TEST_CASE(CommandLineParser_ParsePositional_test) {
	pout << "CommandLineParser_ParsePositional_test..." << endl;

	char arg1[] = "program_name";
	char arg2[] = "pos1";
	char arg3[] = "2";
	char arg4[] = "3.0";
	char arg5[] = "1";
	char arg6[] = "--mogli";
	char arg7[] = "1.0";
	char* argv[] = { arg1, arg2, arg3, arg4, arg5, arg6, arg7, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.setPositionalArgument<string>("arg1", "test");
	c.setPositionalArgument<int>("arg2", "test int");
	c.setPositionalArgument<double>("arg3", "test");
	c.setPositionalArgument<int>("arg4", "test");
	c.setArgument<double>("mogli", "test double", 0.1);
	c.setArgument<double>("def", "test default", 12);
	c.importOptions();

	BOOST_CHECK_EQUAL("pos1", c.getArgument<string>("arg1"));
	BOOST_CHECK_EQUAL(2, c.getArgument<int>("arg2"));
	BOOST_CHECK_EQUAL(3.0, c.getArgument<double>("arg3"));
	BOOST_CHECK_EQUAL(1, c.getArgument<int>("arg4"));
	BOOST_CHECK_EQUAL(1.0, c.getArgument<double>("mogli"));
	BOOST_CHECK_EQUAL(12, c.getArgument<double>("def"));

	pout << "done." << endl;
} /* CommandLineParser_ParsePositional_test */



BOOST_AUTO_TEST_CASE(CommandLineParser_ApplyToConfiguration_test) {
	pout << "CommandLineParser_ApplyToConfiguration_test..." << endl;

	char arg1[] = "program_name";
	char arg2[] = "--standard-lbm";
	char* argv[] = { arg1, arg2, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.setArgument<double>("mogli", "test double", 0.1);
	c.importOptions();

	SolverConfiguration cfg;
	c.applyToSolverConfiguration(cfg);

	BOOST_CHECK_CLOSE(cfg.getCFL(), sqrt(2)*4, 1e-10);
	BOOST_CHECK_EQUAL(cfg.getSedgOrderOfFiniteElement(), 2);

	pout << "done." << endl;
} /* CommandLineParser_ApplyToConfiguration_test */

BOOST_AUTO_TEST_CASE(CommandLineParser_MakePositional_test) {
	pout << "CommandLineParser_MakePositional_test..." << endl;

	char arg1[] = "program_name";
	char arg3[] = "5.0";
	char* argv[] = { arg1, arg3, NULL };
	int argc = sizeof(argv) / sizeof(char*) - 1;

	CommandLineParser c(argc, argv);

	// define arguments
	c.makePositional("cfl");
	c.importOptions();
	BOOST_CHECK_CLOSE(c.getArgument<double>("cfl"), 5.0, 1e-10);

	pout << "done." << endl;
} /* CommandLineParser_MakePositional_test */


BOOST_AUTO_TEST_SUITE_END()

