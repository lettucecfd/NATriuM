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

#include <boost/shared_ptr.hpp>

#include "deal.II/base/logstream.h"

#include "BasicNames.h"

namespace natrium {

enum LogLevel {
	SILENT,
	ERROR,
	WARNING,
	WELCOME,
	BASIC,
	DETAILED,
	DEAL_II_BASIC,
	DEAL_II_DETAILED,
	ALL
};

/**
 * @short this class is responsible for output streams to the command line and log file
 * @note Asn ormally, there is only one ouput stream needed, this class provides the
 * possibily to access the default logging instance from everywhere. The friend functions LOG and LOGGER are used as a
 * wrapper around the private, but static, m_LOGGER. When m_LOGGER does not exists, it is created by LOG or LOGGER.
 */
class Logging: public dealii::LogStream {
	/// friend function; can be used like "LOG(ERROR) << ...;", to write output to the log file and command line
	friend Logging& LOG(LogLevel level);
	/// friend function; is used to access the default logging instance; e.g. LOGGER().setLogFile(...)
	friend Logging& LOGGER();
private:
	/// The Log file has different levels (through sections). This variables stores the current level of the parser.
	size_t m_currentLevel;
	/// the file stream has to be stored as a local variable in order to be not deleted
	shared_ptr<std::ofstream> m_fileStream;
	/// this static object points to the default logstream. It can be accessed through the friend function LOG and LOGGER.
	static shared_ptr<Logging> m_LOGGER;
public:
	/**
	 * @short constructor
	 */
	Logging() {
		m_currentLevel = 0;
		pop();
		setConsoleLevel(BASIC);
		setFileLevel(ALL);
	}

	/// set log level for command line output
	Logging& operator()(LogLevel level) {
		while (m_currentLevel < level) {
			push("");
			m_currentLevel++;
		}
		while (m_currentLevel > level) {
			pop();
			m_currentLevel--;
		}
		return *this;
	}

	/// set log file
	void setLogFile(string logFile) {
		if (has_file()) {
			detach();
		}
		m_fileStream = make_shared<std::ofstream>(logFile);
		attach(*m_fileStream);
	}

	/// detach log file
	void unsetLogFile(string logFile) {
		if (has_file()) {
			detach();
		}
	}

	/// set the level of console output
	LogLevel setConsoleLevel(LogLevel level) {
		return LogLevel(depth_console(level));
	}

	/// set the level of file output
	LogLevel setFileLevel(LogLevel level) {
		return LogLevel(depth_file(level));
	}

};

/**
 *  @short can be used like "LOG(ERROR) << ...;", to write output to the log file and command line
 *  @note accesses the default logstream Logging::m_LOGGER
 */
Logging& LOG(LogLevel level);
/** @short
 *  is used to access the default logging instance m_LOGGER; e.g. LOGGER().setLogFile(...)
 */
Logging& LOGGER();

} /* namespace natrium */

#endif /* LOGGING_H_ */
