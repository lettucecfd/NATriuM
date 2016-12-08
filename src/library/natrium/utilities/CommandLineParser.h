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
 * @short Exception class for CommandLineParser
 */
class HelpMessageStop: public std::exception {
private:
	std::string message;
public:
	HelpMessageStop() {
	}
	~HelpMessageStop() throw () {
	}
	const char *what() const throw () {
		return "";
	}
};

/**
 * @short A class to set NATriuM's solver configuration via the command line.
 *
 * It overwrites default arguments or arguments specified in a configuration file.
 * It also supports user-defined arguments.
 * The usage of this class in a program is demonstrated in step-1. It is used as follows:
 * -# In the beginning of the script, an instance is created: CommandLineParser parse(argc, argv);
 * 	 This has to happen after the call to MPIGuard::getInstance(argc, argv), as the MPIGuard removes the mpirun arguments
 * 	 from argv and argc.
 * -# Then, we can optionally specify a documentation parse.addDocumentationString("some_program_name", "some_program_description")
 * -# Then, the command line arguments are specified
 * 		- parse.setArgument<double>("name", "description", default_value) defines an argument that is interpreted as a double.
 * 		- this argument can be specified in the command line by passing --name <value>
 * 		- the same function exists for double and int
 * 		- NATriuM has a lot of pre-defined arguments, such as 'help'. They may not be overwritten (names have to be unique).
 * 		- parse.setFlag("name", "Description") defines a flag (something like --help)
 * 		- short options can also be specified: parse.setArgument<double>("name,n", "description", default_value)
 * 		- positional arguments are also supported: parse.setPositionalArgument<double>("name", "description")
 * 		  They are defined in the order of their occurence on the command line
 * 		- positional arguments are given to the code like ./program 1.0 0.1
 * -# After all arguments are defined, we call parse.importOptions() to read in the options
 * -# In the code, we can use to the command line arguments by parse.hasArgument() and parse.getArgument<int>("...")
 * -# To apply the arguments to the solver configuration, we have to call parse.applyToSolverConfiguration().
 * 	  This step is recommended to be done after the solver configuration is completely defined, meaning that
 * 	  the command line arguments have higher priority than the configurations that are specified otherwise.
 *
 */
class CommandLineParser {
private:
	int m_argc;
	char** m_argv;
	std::string m_name;
	std::string m_documentation;
	po::options_description m_description;
	po::positional_options_description m_positionalOptions;
	po::variables_map m_varMap;
	bool m_isImported;

	/**
	 * @short Makes sure that the command line arguments are specified in the beginning of the main function.
	 *		  Makes sure that the program stops when the 'help' option is specified.
	 */
	void assert_clean() {
		if (m_isImported) {
			LOG(ERROR)
					<< "Call to CommandLineParser::set... after first CommandLineParser::get... is not allowed."
					<< endl;
			throw CommandLineParserException(
					"Options are already imported. Please define options in the beginning of your program.");
		}
	}

	void assert_imported() {
		if (not m_isImported) {
			LOG(ERROR)
					<< "Call to CommandLineParser::get... is only allowed after CommandLineParser::importOptions() has been called."
					<< endl;
			throw CommandLineParserException(
					"You have to call importOptions before calling a getter function.");
		}
	}

public:

	/// Constructor
	CommandLineParser(int argc, char** argv);

	/// Destructor
	virtual ~CommandLineParser() {

	}

	/**
	 * @short add a global documentation string for the program to the help message
	 */
	void addDocumentationString(std::string program_name, std::string doc) {
		m_name = program_name;
		m_documentation = doc;
	}

	/**
	 * @short import options from argc and argv; Mandatory
	 */
	void importOptions() ;

	/**
	 * @short print the help message
	 */
	friend std::ostream& operator<<(std::ostream& os,
			const CommandLineParser& obj) ;

	/**
	 * @short Define a positional argument of type 'type'
	 */
	template<class type>
	void setPositionalArgument(std::string name, std::string description);

	/**
	 * @short make a non-positional argument positional (e.g. for NATriuM's reserved options)
	 */
	void makePositional(std::string name);

	/**
	 * @short Define a flag (i.e. a command line parameter that comes without value, like --help)
	 */
	void setFlag(std::string name, std::string description);

	/**
	 * @short Define a double parameter (i.e. a command line parameter like --cfl 2.0)
	 */
	template<class type>
	void setArgument(std::string name, std::string description,
			type default_value);

	/**
	 * @short check if the argument is set
	 */
	bool hasArgument(std::string name) {
		assert_imported();
		if (m_varMap.count(name.c_str()))
			return true;
		else
			return false;
	}

	/**
	 * @short get the value of a double argument
	 */
	template <class type>
	type getArgument(std::string name) ;

	/**
	 * @short Overwrite the settings in the solver configuration by the specified command line options
	 */
	void applyToSolverConfiguration(SolverConfiguration& cfg);

}
;

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_COMMANDLINEPARSER_H_ */
