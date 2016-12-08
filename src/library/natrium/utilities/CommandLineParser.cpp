/*
 * CommandLineParser.cpp
 *
 *  Created on: 02.12.2016
 *      Author: akraem3m
 */

#include "CommandLineParser.h"

namespace natrium {

CommandLineParser::CommandLineParser(int argc, char** argv) :
		m_argc(argc), m_argv(argv), m_description("Allowed options"), m_isImported(
				false) {
	m_description.add_options()("help,h", "produce help message");

	// define standard arguments that can be used in all scripts to change the solver configuration
	setFlag("standard-lbm",
			"sets the order to 2 and cfl to 4*sqrt(2), giving standard lbm on regular grids");
	setArgument<int>("order-fe", "order of finite elements, 0 will be ignored", 0);
	setArgument<double>("cfl", "CFL number, 0 will be ignored", 0);
	setFlag("output-on", "switches output on");

}

void CommandLineParser::applyToSolverConfiguration(SolverConfiguration& cfg) {

	// standard lbm on regular grids
	if (hasArgument("standard-lbm")) {
		cfg.setCFL(sqrt(2) * 4);
		cfg.setSedgOrderOfFiniteElement(2);
		LOG(BASIC) << "Scheme set to standard LBM via command line." << endl;
	}

	// finite element order
	assert(hasArgument("order-fe"));
	int order_fe = getArgument<int>("order-fe");
	if (order_fe != 0) {
		cfg.setSedgOrderOfFiniteElement(order_fe);
		LOG(BASIC) << "Order set to " << order_fe << " via command line."
				<< endl;
	}

	// cfl number
	assert(hasArgument("cfl"));
	double cfl = getArgument<double>("cfl");
	if (cfl != 0) {
		cfg.setCFL(cfl);
		LOG(BASIC) << "CFL set to " << cfl << " via command line." << endl;
	}

	// switch output on
	if (hasArgument("output-on")) {
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Output switched on via command line." << endl;
	}
}

template<class type>
void CommandLineParser::setPositionalArgument(std::string name,
		std::string description) {
	assert_clean();
	try {
		m_description.add_options()(name.c_str(), po::value<type>(),
				description.c_str());
		m_positionalOptions.add(name.c_str(), 1);
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the positional argument ";
		msg << name;
		msg << "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "The following arguments were set previously:" << endl;
		msg << m_description << endl;
		throw CommandLineParserException(msg.str());
	}
}

void CommandLineParser::makePositional(std::string name) {
	assert_clean();
	/*if (not hasArgument(name)){
		std::stringstream msg;
		msg << "A problem occured when making ";
		msg << name;
		msg << " a positional argument." << endl;
		msg << endl;
		msg << "The argument '" << name << "' was not specified before the call to makePositional()." << endl;
		msg << m_description << endl;
		throw CommandLineParserException(msg.str());

	}*/
	try {
		m_positionalOptions.add(name.c_str(), 1);
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the positional argument ";
		msg << name;
		msg << "Make sure that the name of the argument was defined (but not as a positional argument) before." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "Help message:" << endl;
		msg << *this << endl;
		throw CommandLineParserException(msg.str());
	}
}

//explicit instantiation
template
void CommandLineParser::setPositionalArgument<int>(std::string name,
		std::string description);
template
void CommandLineParser::setPositionalArgument<double>(std::string name,
		std::string description);
template
void CommandLineParser::setPositionalArgument<string>(std::string name,
		std::string description);

template<class type>
void CommandLineParser::setArgument(std::string name, std::string description,
		type default_value) {
	assert_clean();
	try {
		type opt;
		m_description.add_options()(name.c_str(),
				po::value<type>(&opt)->default_value(default_value),
				description.c_str());
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the argument ";
		msg << name;
		msg << "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "The following arguments were set previously:" << endl;
		msg << m_description << endl;
		throw CommandLineParserException(msg.str());
	}
}

template void CommandLineParser::setArgument<int>(std::string name,
		std::string description, int default_value);
template void CommandLineParser::setArgument<double>(std::string name,
		std::string description, double default_value);
template void CommandLineParser::setArgument<string>(std::string name,
		std::string description, string default_value);

template<class type>
type CommandLineParser::getArgument(std::string name) {
	assert_imported();
	if (m_varMap.count(name)) {
		try {
			return m_varMap[name].as<type>();
		} catch (po::error& e) {
			std::stringstream msg;
			msg << "A problem occured when reading the argument ";
			msg << name;
			msg << " (of type double.)"
					<< "Make sure that the argument is properly defined.";
			msg << endl;
			msg << "Original error message:" << e.what();
			throw CommandLineParserException(msg.str());
		}
	} else {
		throw CommandLineParserException(
				"Did not find command line argument. Pass --help to see the help message.");
	}
}
template double CommandLineParser::getArgument<double>(std::string name);
template int CommandLineParser::getArgument<int>(std::string name);
template string CommandLineParser::getArgument<string>(std::string name);

void CommandLineParser::setFlag(std::string name, std::string description) {
	assert_clean();
	try {
		m_description.add_options()(name.c_str(), description.c_str());
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the flag ";
		msg << name;
		msg << " (of type flag.)"
				<< "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		throw CommandLineParserException(msg.str());
	}
}

void CommandLineParser::importOptions() {

	m_isImported = true;

	// positional arguments
	po::store(
			po::command_line_parser(m_argc, m_argv).options(m_description).positional(
					m_positionalOptions).run(), m_varMap);
	po::notify(m_varMap);

	// print out help message and stop program
	if (m_varMap.count("help")) {
		pout << *this << endl;
		throw HelpMessageStop();
	}
}

std::ostream& operator<<(std::ostream& os,
		const CommandLineParser& obj) {
	os << endl;
	os << "Help message for NATriuM executable " << obj.m_name << endl;
	os << "------------------------------------------------------" << endl;
	os << "Usage: \t" << obj.m_name << "  ";
	for (size_t i = 0; i < obj.m_positionalOptions.max_total_count(); i++) {
		os << "<" << obj.m_positionalOptions.name_for_position(i) << "> ";
	}
	os << endl;
	os << endl;
	os << "Documentation:" << endl;
	os << obj.m_documentation << endl;
	os << "NATriuM is a library for off-lattice Boltzmann simulations."
			<< endl;
	os << endl;
	os << obj.m_description << endl;
	return os;
}


} /* namespace natrium */

