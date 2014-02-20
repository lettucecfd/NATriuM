/**
 * @file Logging.h
 * @short Definition of logging output streams
 * @date 19.02.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Logging.h"

#include <fstream>
#include <sstream>

#include "boost/make_shared.hpp"

namespace natrium {

// use only shared pointers: otherwise this stuff is never removed again

boost::shared_ptr<std::ofstream> fileToNull = boost::make_shared<std::ofstream>(
		"../results/test/natrium.log");

boost::shared_ptr<std::stringstream> logToNull = boost::make_shared<
		std::stringstream>();

// Create Tee objects, which allow output to file and cout, simultaneously
Tee fullTee(
		*logToNull, *fileToNull);
Tee basicTee(
		std::cout, *fileToNull);
Tee errorTee(
		std::cerr, *fileToNull);

// Create static Stream objects
boost::shared_ptr<TeeStream> Logging::FULL= boost::make_shared<TeeStream>(fullTee);
boost::shared_ptr<TeeStream> Logging::BASIC= boost::make_shared<TeeStream>(basicTee);
boost::shared_ptr<TeeStream> Logging::ERROR= boost::make_shared<TeeStream>(errorTee);

}
/* namespace natrium */
