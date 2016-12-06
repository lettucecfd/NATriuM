/**
 * @short a module to read arguments from the command line
 */

#ifndef LIBRARY_NATRIUM_COMMANDLINEPARSER_H_
#define LIBRARY_NATRIUM_COMMANDLINEPARSER_H_

#include "boost/program_options.hpp"
#include "BasicNames.h"
#include "../solver/SolverConfiguration.h"

namespace po = boost::program_options;
// for instructions, check out: http://www.boost.org/doc/libs/1_41_0/doc/html/program_options/tutorial.html

namespace natrium {



/**
 * @short Exception class for CommandLineParser
 */
class CommandLineParserException: public NATriuMException {
private:
	std::string message;
public:
	CommandLineParserException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	CommandLineParserException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~CommandLineParserException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short A class to set NATriuM's solver configuration via the command line,
 * 			overwriting default arguments or arguments specified in a file.
 * 			It also supports user arguments.
 */
class CommandLineParser {
private:
	int m_argc;
	char** m_argv;
	po::options_description m_description;
	po::variables_map m_varMap;
	bool m_isImported;

	/**
	 * @short import options from argc and argv; is private and implemented in a way that the user can just call the getter functions
	 * 			and the CommandLineParser does the rest automatically
	 */
	void import_options() {
		po::store(po::parse_command_line(m_argc, m_argv, m_description),
				m_varMap);
		po::notify(m_varMap);
	}

public:

	/// Constructor
	CommandLineParser(int argc, char** argv);

	void setDoubleArgument(std::string name, double default_value,
			size_t position, std::string description) {

	}

	void setIntArgument(std::string name, int default_value, size_t position,
			std::string description) {
		m_description.add_options()
				(name.c_str(), po::value<int>(), description.c_str());

	}

	void setStringArgument(std::string name, string default_value,
			size_t position, std::string description) {

	}

	double getDoubleArgument(std::string name) {
		if (not m_isImported) {
			import_options();
			m_isImported = true;
		}
		if (m_varMap.count(name)){
			return m_varMap[name].as<double>();
		} else {
			throw CommandLineParserException("Did not find command line argument");
		}

	}

	int getIntArgument(std::string name) {
		if (not m_isImported) {
			import_options();
			m_isImported = true;
		}
		if (m_varMap.count(name)){
			return m_varMap[name].as<int>();
		} else {
			throw CommandLineParserException("Did not find command line argument");
		}
	}

	std::string getStringArgument(std::string name) {
		if (not m_isImported) {
			import_options();
			m_isImported = true;
		}
		if (m_varMap.count(name)){
			return m_varMap[name].as<std::string>();
		} else {
			throw CommandLineParserException("Did not find command line argument");
		}
	}

	void applyToSolverConfiguration(SolverConfiguration& cfg);

}
;

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_COMMANDLINEPARSER_H_ */
