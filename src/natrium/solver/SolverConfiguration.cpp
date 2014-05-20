/**
 * @file SolverConfiguration.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SolverConfiguration.h"

namespace natrium {

SolverConfiguration::SolverConfiguration() {
	// Declare structure of the parameter file
	enter_subsection("General");
	{
		declare_entry("Time step size", "0.2",
				dealii::Patterns::Double(1e-10),
				"Size of the (initial) time step.");
		declare_entry("Stencil", "D2Q9",
				dealii::Patterns::Selection("D2Q9"),
				"The discrete velocity stencil. The number behind D denotes the dimension (2 or 3). The number behind Q denotes the number of particle directions in the discrete velocity model.");
		declare_entry("Stencil scaling", "1.0",
				dealii::Patterns::Double(1e-10),
				"The scaling of the discrete velocities. Whereas in the standard LBM the magnitude of the particle velocities is set to 1.0 due to the uniform mesh grid, the SEDG-LBM features scaled particle velocities. As the scaling factor is proportional to the speed of sound, it strongly impacts the relaxation time.");
	}
	leave_subsection();
	enter_subsection("Advection");
	{
		declare_entry("Advection scheme", "SEDG",
				dealii::Patterns::Selection("SEDG"),
				"The algorithm which is used for the advection (=streaming) step. While the LBM on a uniform mesh facilitates streaming towards a simple index shift, non-uniform meshes need a more sophisticated advection scheme.");
		declare_entry("Time integrator", "Runge-Kutta 5-stage",
				dealii::Patterns::Selection("Runge-Kutta 5-stage"),
				"The algorithm which is used for the time integration of the discretizted advection (=streaming) equation. A time integrator is required, when the advection scheme is based upon some Finite Element/Difference/Volume or discontinuous Galerkin scheme.");

		enter_subsection("SEDG");
		{
			declare_entry("Order of finite element", "4",
					dealii::Patterns::Integer(2),
					"The degree of the polynomial shape functions used by the SEDG scheme.");
			declare_entry("Flux type", "Lax-Friedrichs",
					dealii::Patterns::Selection("Lax-Friedrichs|Central"),
					"The flux connects the shape functions between neighboring elements. It is strongly recommended to use the Lax-Friedrichs scheme which is a forward-discretization along characteristics, rather than a central flux.");
		}
		leave_subsection();
	}
	leave_subsection();
	enter_subsection("Collision");
	{
		declare_entry("Collision scheme",
				"BGK with transformed distribution functions",
				dealii::Patterns::Selection(
						"BGK with transformed distribution functions"),
				"The collision step models velocity changes due to particle collisions (local at each node) by a relaxation towards thermodynamic equilibrium. There are several approaches, e.g. the single-relaxation time Bhatnagar-Groß-Krook. Using transformed particle distribution functions enhances the accuracy of the LBM.");
		declare_entry("Collision on boundary nodes", "true",
				dealii::Patterns::Bool(),
				"States whether the collision step is to be done on all nodes or only on internal nodes. E.g. the standard bounce back scheme is of 2nd order, when collisions take place at boundary nodes, and of 1st order, if not.");

	}
	leave_subsection();

	enter_subsection("Initialization");
	{
		declare_entry("Restart at last checkpoint?", "false",
				dealii::Patterns::Bool(),
				"The solver can be restarted at the last stored checkpoint, in case that an old run had been aborted at some point of time.");
		declare_entry("Initialization scheme", "Equilibrium",
				dealii::Patterns::Selection("Equilibrium|Iterative"),
				"The initial particle distribution functions are normally assumed to be in local equilibrium. A more stable (and costly) scheme is to do some streaming steps on the density field but not on the velocity field, before starting the actual simulations (see e.g. the Book of Guo and Shu).");
		enter_subsection("Iterative initialization stop condition");
		{
			declare_entry("Residual", "1e-6",
					dealii::Patterns::Double(1e-25),
					"The iterative initialization stops, when the density increment is smaller than the residual, i.e. the iteration has converged.");
			declare_entry("Number of iterations", "2000",
					dealii::Patterns::Integer(1),
					"The iterative initialization stops at the latest after a specific number of iterations.");

		}
		leave_subsection();
	}
	leave_subsection();

	enter_subsection("Stop condition");
	{
		declare_entry("Number of time steps", "1000000",
				dealii::Patterns::Integer(1),
				"The maximum number of time steps.");

	}
	leave_subsection();

	enter_subsection("Output");
	{
		declare_entry("Switch output off?", "false",
				dealii::Patterns::Bool(), "Switch output off, completely.");
		declare_entry("Output directory", "/tmp/NATriuM",
				dealii::Patterns::DirectoryName(),
				"The name of the directory to which the output is written.");
		declare_entry("Output checkpoint interval", "1000",
				dealii::Patterns::Integer(1),
				"Write out checkpoint files every ... step.");
		declare_entry("Output solution interval", "1000",
				dealii::Patterns::Integer(1),
				"Write out solution every ... step.");
		declare_entry("Command line verbosity", "Basic",
				dealii::Patterns::Selection("Error|Basic|Full"),
				"The amount of command line output.");
		declare_entry("Write a log file?", "true", dealii::Patterns::Bool(),
				"Specifies if log is written to a file.");
	}
	leave_subsection();


}

SolverConfiguration::SolverConfiguration(const std::string& XMLfilename){
	SolverConfiguration();
	readFromXMLFile(XMLfilename);
}


/**
 * @short wrapper function for ParameterHandler::read_input; directing cerr into a C++-Exception
 **/
void SolverConfiguration::readFromTextFile(const std::string & filename, const bool optional, const bool write_stripped_file) {
	// redirect cerr to a string buffer
	std::stringstream buffer;
	std::streambuf * old = cerr.rdbuf(buffer.rdbuf());
	bool isEverythingOK = ParameterHandler::read_input(filename, optional,
			write_stripped_file);
	std::string errorMessage = buffer.str();
	// reset cerr to console output
	cerr.rdbuf(old);

	if (not isEverythingOK) {
		throw ConfigurationException(errorMessage);
	}
} /* readFromTextFile */

/**
 * @short wrapper function for ParameterHandler::read_input_from_xml; directing cerr into a C++-Exception
 **/
void SolverConfiguration::readFromXMLFile(const std::string & filename) {
	// redirect cerr to a string buffer
	std::stringstream buffer;
	std::streambuf * old = cerr.rdbuf(buffer.rdbuf());

	// create stream
	std::ifstream xmlFile(filename);
	bool isEverythingOK = ParameterHandler::read_input_from_xml(xmlFile);
	std::string errorMessage = buffer.str();
	// reset cerr to console output
	cerr.rdbuf(old);

	if (not isEverythingOK) {
		throw ConfigurationException(errorMessage);
	}
} /* readFromXMLFile */


} /* namespace natrium */
