/**
 * @file Logging_test.cpp
 * @short
 * @date 23.05.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/Logging.h"

#include <time.h>

#include "boost/test/unit_test.hpp"


//#define TEST_LOGGING
#ifdef TEST_LOGGING

using namespace natrium;

BOOST_AUTO_TEST_SUITE(Logging_test)

BOOST_AUTO_TEST_CASE(Logging_Levels_of_Parent_Class_test){
	cout << "Logging_Levels_of_Parent_Class_test..." << endl;
	Logging a;
	a.pop();
	BOOST_CHECK(not a.has_file());
	a.attach(cout, true);
	// cout functions as a file
	BOOST_CHECK(a.has_file());
	a << "This message should go to cout and cerr." << endl;
	// Test different output on different levels
	a.depth_console(1);
	a.depth_file(2);
	a.push("");
	a << "This message should go to cout and cerr." << endl;
	a.push("");
	a << "This message should go to cout only." << endl;
	a.pop();
	a.pop();
	// test time output
	sleep(1);
	// For some reason swapping the next 2 lines is erroneous
	a.log_cerr();
	a.detach();
	cout << "guppy" << endl;
	cout << "guppY" << endl;
	a.log_execution_time(true);
	a.depth_console(5);
	a << "This is the time." << endl;
	a.log_time_differences(true);
	sleep(1);
	a << "Time passed since last log" << endl;
	a.log_execution_time(false);
	a.log_time_differences(false);
	a << "Printing a double: " << 0.0001 << endl;
	cout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(Logging_DifferentLevelsWrapper_test){
	cout << "Logging_DifferentLevelsWrapper_test..." << endl;

	Logging a;

	a.depth_console(ERROR);
	a(ERROR) << "This is an Error" << endl;
	a(WARNING) << "The warning is not displayed" << endl;

	cout << "done." << endl;
}

BOOST_AUTO_TEST_CASE(Logging_WithGlobalObjectLOG_test){
	LOG(ERROR) << "This should go only to stderr (global object)" << endl;
	LOGGER().setConsoleLevel(ERROR);
	LOG(WARNING) << "The warning should not appear." << endl;
	LOGGER().setConsoleLevel(WARNING);
	LOG(WARNING) << "Now, the warning should appear." << endl;
	LOGGER().setLogFile("/tmp/natrium.log");
	LOG(ERROR) << "This should go to the file and console." << endl;

}


BOOST_AUTO_TEST_SUITE_END()

#endif




