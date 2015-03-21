/**
 * @file SolverConfiguration.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SolverConfiguration.h"

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

namespace natrium {

SolverConfiguration::SolverConfiguration() {
	// Declare structure of the parameter file
	enter_subsection("General");
	{
		declare_entry("Time step size", "0.2", dealii::Patterns::Double(1e-10),
				"Size of the (initial) time step.");
		declare_entry("Stencil", "D2Q9", dealii::Patterns::Selection("D2Q9"),
				"The discrete velocity stencil. The number behind D denotes the dimension (2 or 3). The number behind Q denotes the number of particle directions in the discrete velocity model.");
		declare_entry("Stencil scaling", "1.0", dealii::Patterns::Double(1e-10),
				"The scaling of the discrete velocities. Whereas in the standard LBM the magnitude of the particle velocities is set to 1.0 due to the uniform mesh grid, the SEDG-LBM features scaled particle velocities. As the scaling factor is proportional to the speed of sound, it strongly impacts the relaxation time.");
		declare_entry("Has analytic solution?", "false",
				dealii::Patterns::Bool(),
				"Indicates whether the given flow problem is inherited of BenchmarkProblem, which means that it has an analytic solution with which the numerical solution can be compared.");
	}
	leave_subsection();
	enter_subsection("Advection");
	{
		declare_entry("Advection scheme", "SEDG",
				dealii::Patterns::Selection("SEDG"),
				"The algorithm which is used for the advection (=streaming) step. While the LBM on a uniform mesh facilitates streaming towards a simple index shift, non-uniform meshes need a more sophisticated advection scheme.");
		declare_entry("Time integrator", "Runge-Kutta 5-stage",
				dealii::Patterns::Selection(
						"Runge-Kutta 5-stage|Theta method|Exponential|Other"),
				"The algorithm which is used for the time integration of the discretizted advection (=streaming) equation. "
						"A time integrator is required, when the advection scheme is based upon some Finite Element/Difference/Volume "
						"or discontinuous Galerkin scheme. Other refers to  deal.II's built-in integrators which are accessible through "
						"the section DealIntegrator.");
		enter_subsection("SEDG");
		{
			declare_entry("Order of finite element", "1",
					dealii::Patterns::Integer(1),
					"The degree of the polynomial shape functions used by the SEDG scheme.");
			declare_entry("Flux type", "Lax-Friedrichs",
					dealii::Patterns::Selection("Lax-Friedrichs|Central"),
					"The flux connects the shape functions between neighboring elements. It is strongly recommended to use the Lax-Friedrichs scheme which is a forward-discretization along characteristics, rather than a central flux.");
		}
		leave_subsection();
		enter_subsection("Theta method");
		{
			declare_entry("Theta", "0.5", dealii::Patterns::Double(0, 1),
					"Theta=0: Explicit Euler; Theta=0.5: Crank-Nicolson; Theta=1.0: Implicit Euler.");
		}
		leave_subsection();
		enter_subsection("Deal.II integrator");
		{
			declare_entry("Runge Kutta scheme", "None",
					dealii::Patterns::Selection(
							"Forward Euler|RK 3rd order|RK Classic 4th order|Backward Euler|"
									"Implicit midpoint|Crank-Nicoloson|SDIRK 2 stages|Heun-Euler|Bogacki-Shampine|Dopri|Fehlberg|Cash-Karp|None"),
					"Deal.ii built-in time integrators. They are accessed by chosing time integrator = other.");
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
			declare_entry("Residual", "1e-6", dealii::Patterns::Double(1e-25),
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
		declare_entry("Switch output off?", "false", dealii::Patterns::Bool(),
				"Switch output off, completely.");
		declare_entry("User interaction?", "true", dealii::Patterns::Bool(),
				"Indicates, whether to ask the user in case that configurations are unclear or old files are possibly overwritten.");
		declare_entry("Output directory", "/tmp/NATriuM",
				dealii::Patterns::DirectoryName(),
				"The name of the directory to which the output is written.");
		declare_entry("Output checkpoint interval", "1000",
				dealii::Patterns::Integer(1),
				"Write out checkpoint files every ... step.");
		declare_entry("Output table interval", "1000",
				dealii::Patterns::Integer(1),
				"Write out a line to the result table every ... step.");
		declare_entry("Output solution interval", "1000",
				dealii::Patterns::Integer(1),
				"Write out solution every ... step.");
		declare_entry("Command line verbosity", "5",
				dealii::Patterns::Integer(0, 8),
				"The amount of command line output.");
		declare_entry("Write a log file?", "true", dealii::Patterns::Bool(),
				"Specifies if log is written to a file.");
	}
	leave_subsection();

}

SolverConfiguration::SolverConfiguration(const std::string& XMLfilename) {
	SolverConfiguration();
	readFromXMLFile(XMLfilename);
}

/**
 * @short wrapper function for ParameterHandler::read_input; directing cerr into a C++-Exception
 **/
void SolverConfiguration::readFromTextFile(const std::string & filename,
		const bool optional, const bool write_stripped_file) {
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

void SolverConfiguration::prepareOutputDirectory() {
	/// If not exists, try to create output directory
	//  ((Using boost::filesystem provides a cross-platform solution))
	boost::filesystem::path outputDir(getOutputDirectory());
	boost::filesystem::path parentDir(outputDir.branch_path());
	if (not boost::filesystem::is_directory(parentDir)) {
		std::stringstream msg;
		msg << "You want to put your output directory into "
				<< parentDir.string()
				<< ", but this parent directory does not even exist.";
		throw ConfigurationException(msg.str());
	}
	// Make output directory
	try {
		//create_directory throws basic_filesystem_error<Path>, if fail other than that the directory already existed
		//returns false, if directory already existed
		bool newlyCreated = boost::filesystem::create_directory(outputDir);
		if (newlyCreated) {
			return;
		}
	} catch (std::exception& e) {
		std::stringstream msg;
		msg << "You want to put your output directory into "
				<< parentDir.string()
				<< ", but you seem to have no writing permissions.";
		throw ConfigurationException(msg.str());
	}
	// Postcondition: directory exists
	// Check writing permissions in directory
	try {
		/// try to create a single file
		std::ofstream filestream;
		filestream.open((outputDir / "testtatata.txt").string().c_str());
		filestream << " ";
		filestream.close();
		boost::filesystem::remove(
				(outputDir / "testtatata.txt").string().c_str());
		/// try to open all files
		boost::filesystem::directory_iterator it(outputDir), eod;
		BOOST_FOREACH(boost::filesystem::path const &p, std::make_pair(it, eod)) {
			if (not boost::filesystem::is_directory(p)) {
				std::fstream filestream;
				filestream.open(p.string().c_str(),
						std::fstream::app | std::fstream::out);
				// throw exception if file is not opened
				if (not filestream.is_open()) {
					throw std::exception();
				}
			}
		}
	} catch (std::exception& e) {
		std::stringstream msg;
		msg
				<< "You don't have writing access to the files which are already existing in your Output directory "
				<< outputDir.string();
		throw ConfigurationException(msg.str());
	}
	// check if something is possibly going to be overwritten
	clock_t begin = clock();
	if ((not isRestartAtLastCheckpoint())
			and (not boost::filesystem::is_empty(outputDir))) {
		if (isUserInteraction()) {
			// Request user input
			cout
					<< "'Restart at last checkpoint' is disabled, but Output directory is not empty. The simulation might overwrite old data. Do you really want to continue?"
					<< endl;
			size_t yes1_or_no2 = 0; // = 1 for yes; = 2 for no
			string input = "";
			for (size_t i = 0; true; i++) {
				cout << "Please enter 'y' or 'n':" << endl;
				getline(std::cin, input);
				// check for yes
				if ("y" == input) {
					yes1_or_no2 = 1;
					break;
					// check for no
				} else if ("n" == input) {
					yes1_or_no2 = 2;
					break;
				}
				// check for too many tries
				if (i > 5) {
					break;
				}
				// check for timeout
				clock_t end = clock();
				double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				if (elapsed_secs > 30) {
					break;
				}
				cout << "Your input was not understood. ";
			}
			// no sound input
			if (0 == yes1_or_no2) {
				throw ConfigurationException(
						"Requested user input, but did not get meaningful answer.");
			} else if (2 == yes1_or_no2) {
				throw ConfigurationException(
						"Execution stopped due to user's intervention.");
			}
		} else {
			LOG(BASIC) << "Starting NATriuM..." << endl;
			LOG(WARNING)
					<< "Simulation might overwrite old data in output file."
					<< endl;
		}
	}
}

void SolverConfiguration::isConsistent() {
	// check consistency of time integrator setting
	if (OTHER == this->getTimeIntegrator()) {
		if (NONE == this->getDealIntegrator()) {
			std::stringstream msg1;
			msg1
					<< "If you set time integrator to 'other' you will need to specify the Deal integrator. Found none.";
			LOG(ERROR) << msg1.str().c_str() << endl;
			throw ConfigurationException(msg1.str());
		}
	} else if (NONE != this->getDealIntegrator()) {
		// Warn if the Deal.II integrator setting will not be applied
		LOG(WARNING)
				<< "Did not understand setting of Deal.II integrator. If you want to use the Deal.II "
						"time integration schemes, you will have to set Time integrator to 'OTHER'." << endl;
	}
}

} /* namespace natrium */
