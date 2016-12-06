/*
 * CommandLineParser.cpp
 *
 *  Created on: 02.12.2016
 *      Author: akraem3m
 */

#include "CommandLineParser.h"

namespace natrium {


CommandLineParser::CommandLineParser(int argc, char** argv) :
		//	"This executable is compiled with NATriuM, a package for off-lattice Boltzmann simulations.
		m_argc(argc), m_argv(argv), m_description("Allowed options"), m_isImported(false) {
	m_description.add_options()("help", "produce help message");
}


void CommandLineParser::applyToSolverConfiguration(SolverConfiguration& cfg){

}

} /* namespace natrium */


