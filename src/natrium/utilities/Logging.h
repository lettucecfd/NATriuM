/**
 * @file Logging.h
 * @short Definition of logging output streams
 * @date 18.02.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef LOGGING_H_
#define LOGGING_H_


#include <fstream>
#include <iostream>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/shared_ptr.hpp>

#include "BasicNames.h"

namespace natrium {

typedef boost::iostreams::tee_device<std::ostream, std::ofstream>  Tee;
typedef boost::iostreams::stream< Tee > TeeStream;

enum LogLevel {
	SILENT, ERROR, WARNING, WELCOME, BASIC, DETAILED, DEAL_II_BASIC, DEAL_II_DETAILED, ALL
};


/**
 * @short this class is responsible for output streams to the command line and log file
 */
class Logging {

public:
	/// set log level for command line output

	/// print
	static void print(LogLevel level, string msg);
	/// Full (complete) log; stream for detailed information
	static boost::shared_ptr<TeeStream> FULL;
	/// Stream for basic information
	static boost::shared_ptr<TeeStream> BASIC;
	/// Stream for error messages
	static boost::shared_ptr<TeeStream> ERROR;
};

} /* namespace natrium */

#endif /* LOGGING_H_ */
