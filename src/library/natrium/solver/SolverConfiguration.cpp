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
		declare_entry("CFL", "0.4", dealii::Patterns::Double(1e-10),
				"CFL number. Determines the size of the (initial) time step. The CFL number is defined as stencil_scaling/(dx*(p+1)^2).");
        declare_entry("Prandtl", "1.0", dealii::Patterns::Double(1e-10),
                      "Prandtl number. Determines ratio of viscosity and heat conduction. Set to Pr=0.7 for air. Default Pr=1.0 (BGK)");
        declare_entry("Gamma", "1.4", dealii::Patterns::Double(1e-10),
                      "Heat capacity ratio . Set to gamma=1.4 for air (Default)");
        declare_entry("Prandtl number set", "false",
                      dealii::Patterns::Bool(),
                      "Indicates if Prandtl number deviates from 1.0");
		declare_entry("Stencil", "D2Q9",
				dealii::Patterns::Selection("D2Q9|D3Q13|D3Q19|D3Q15|D3Q21|D3Q27|RD3Q27|D2Q19V|D2Q19H|D2Q25H|D3Q45|D3Q77|D3V27"),
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
				dealii::Patterns::Selection("SEDG|Semi-Lagrangian"),
				"The algorithm which is used for the advection (=streaming) step. While the LBM on a uniform mesh facilitates streaming towards a simple index shift, non-uniform meshes need a more sophisticated advection scheme. SEDG stands for spectral element discontinuous Galerkin. Note that the Semi-Lagrangian streaming does not require a time integrator.");

		declare_entry("Support points", "Gauss-Lobatto",
				dealii::Patterns::Selection(
						"Gauss-Lobatto|Gauss-Lobatto-Chebyshev|Gauss-Chebyshev|Equidistant"),
				"The support points of the finite elements. For the SEDG streaming, "
						"this has to be set to 'Gauss-Lobatto'. For the semi-Lagrangian "
						"streaming, this choice affects the stability for p>2. "
						"We recommend to use Gauss-Lobatto-Chebyshev.");

		declare_entry("Quadrature", "Gauss-Lobatto",
				dealii::Patterns::Selection("Gauss-Lobatto|Gauss"),
				"The quadrature rule that is used for integration of the finite elements. "
						"For the SEDG streaming, this has to be set to 'Gauss-Lobatto'.");

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
		enter_subsection("Deal.II linear solver");
		{
			declare_entry("Linear solver", "Bicgstab",
					dealii::Patterns::Selection(
							"Bicgstab|Cg|Fgmres|Gmres|Minres|Qmrs|Relaxation|Richardson"),
					"Deal.ii built-in linear solvers.");
		}
		leave_subsection();
		enter_subsection("Embedded Parameters");
		{
			declare_entry("Coarsen parameter", "1.2",
					dealii::Patterns::Double(1, 1000),
					"Parameter for the embedded deal.II methods. This parameter is the factor (>1) by which the time step is multiplied when the time stepping can be coarsen.");
			declare_entry("Refinement parameter", "0.8",
					dealii::Patterns::Double(0, 1),
					"Parameter for the embedded deal.II methods. This parameter is the factor (<1) by which the time step is multiplied when the time stepping must be refined.");
			declare_entry("Minimum CFL", "10", dealii::Patterns::Double(),
					"Parameter for the embedded deal.II methods. Smallest CFL allowed.");
			declare_entry("Maximum CFL", "0.05", dealii::Patterns::Double(),
					"Parameter for the embedded deal.II methods. Largest CFL allowed.");
			declare_entry("Refinement tolerance", "1e-8",
					dealii::Patterns::Double(),
					"Parameter for the embedded deal.II methods. Refinement tolerance: if the error estimate is larger than refine_tol, the time step is refined.");
			declare_entry("Coarsen tolerance", "1e-12",
					dealii::Patterns::Double(),
					"Parameter for the embedded deal.II methods. Coarsening tolerance: if the error estimate is smaller than coarse_tol, the time step is coarsen.");
		}
		leave_subsection();
	}
	leave_subsection();
	enter_subsection("Collision");
	{
		declare_entry("Collision scheme", "BGK standard",
				dealii::Patterns::Selection(
						"BGK standard|BGK steady state|BGK standard transformed|BGK multiphase|BGK incompressible|BGK regularized|MRT standard|MRT entropic|Entropic stabilized|KBC standard|KBC central|BGK multi am4|BGK multi bdf2"),
				"The collision step models velocity changes due to particle collisions (local at each node) by a relaxation towards "
						"thermodynamic equilibrium. There are several approaches, e.g. the single-relaxation time Bhatnagar-Groß-Krook (BGK) model. "
						"The standard");

		declare_entry("MRT basis", "Dellar D2Q9",
				dealii::Patterns::Selection(
						"Dellar D2Q9|Lallemand D2Q9|DHumieres D3Q19"),
				"Moment basis for MRT collisions");

		declare_entry("MRT relaxation times", "Full",
				dealii::Patterns::Selection(
						"Full|Dellar D2Q9 Only N|DHumieres Paper"),
				"Relaxation scheme (choice of the higher-order relaxation parameters). Default: full relaxation to equilibrium");

		declare_entry("Equilibrium scheme", "BGK equilibrium",
				dealii::Patterns::Selection("BGK equilibrium|Quartic equilibrium|Incompressible equilibrium|Steady-state equilibrium|Entropic equilibrium"),
				"Defines the equilibrium for the collision.");

		enter_subsection("BGK parameters");
		{
			declare_entry("Steady state gamma", "1",
					dealii::Patterns::Double(0, 1 + 1e-50),
					"The parameter of the steady state preconditioner. For gamma = 1, the scheme is equivalent to the standard BGK"
							"For gamma -> 0, the convergence to steady states is speed up and the effective Mach number is lowered, which"
							"gives nearly incompressible results.");
			declare_entry("Pseudopotential type", "ShanChen",
					dealii::Patterns::Selection(
							"ShanChen|Sukop|CarnahanStarling"),
					"The functional form of the pseudopotential in multiphase simulations.");
			declare_entry("Pseudopotential G", "-5",
					dealii::Patterns::Double(-1e10, 1e10),
					"The parameter that describes the interaction strength in the Pseudopotential multiphase model.");
			declare_entry("Pseudopotential T", "0.0848997582",
					dealii::Patterns::Double(-1e10, 1e10),
					"The parameter that describes the temperature in the Pseudopotential multiphase model (only for Carnahan-Starling EOS).");
			declare_entry("Forcing scheme", "No Forcing",
					dealii::Patterns::Selection(
							"No Forcing|Shifting Velocity|Exact Difference|Guo"),
					"The way to incorporate the force into the LBGK equation.");

		}
		leave_subsection();

	}
	leave_subsection();

	enter_subsection("Filtering");
	{
		declare_entry("Regularization", "No Regularization",
				dealii::Patterns::Selection(
						"No Regularization|Pseudo-entropy maximization|Zero high-order moments|Entropy maximization|Pseudo-entropy maximization with e"),
				"Regularization of high-order moments.");
		declare_entry("Apply vmult limiter?", "false", dealii::Patterns::Bool(),
				"Limiting suppresses oscillations in the semi-Lagrangian streaming "
						"(makes sense when the order of finite elements is > 2).");
		declare_entry("Apply filtering?", "false", dealii::Patterns::Bool());
		declare_entry("Filtering scheme", "Exponential",
				dealii::Patterns::Selection("Exponential|New"),
				"A filter that dampens the high-frequent oscillations from the distribution functions.");
		enter_subsection("Filter parameters");
		{
			declare_entry("Filter interval", "1",
					dealii::Patterns::Integer(1, 10000000),
					"The filter is applied each ... time step.");
			declare_entry("Exponential alpha", "10.0",
					dealii::Patterns::Double(0, 1e10),
					"The exponential filter is defined exp(-alpha * ((poly_degree + 1 -Nc) / (max_poly_degree + 1 -Nc)) ^ s");
			declare_entry("Exponential s", "20.0",
					dealii::Patterns::Double(0, 1e10),
					"The exponential filter is defined exp(-alpha * ((poly_degree + 1 -Nc) / (max_poly_degree + 1 -Nc)) ^ s");
			declare_entry("Exponential Nc", "1",
					dealii::Patterns::Double(0, 50),
					"First polynomial degree that is filtered in the exponential filter, Nc.");
			declare_entry("Degree by component sums", "false",
					dealii::Patterns::Bool(),
					"If true, the degree of the polynomial is calculated as the sum of the one-dimensional degrees, "
							"which creates a much more dissipative filter."
							"Otherwise (default) the degree is calculated as the maximum degree of the one-dimensional shape functions.");
		}
		leave_subsection();
	}
	leave_subsection();

	enter_subsection("Initialization");
	{
		declare_entry("Restart at iteration", "0", dealii::Patterns::Integer(0),
				"The solver can be restarted at a stored checkpoint, in case that an old run had been aborted at some point of time."
						"The iteration at which the solver is to be restarted."
						" You have to make sure that a checkpoint file corresponding to that iteration exists.");
		declare_entry("Initialization scheme", "Equilibrium",
				dealii::Patterns::Selection("Equilibrium|Iterative|CompressibleIterative|Gradients"),
				"The initial particle distribution functions are normally assumed to be in local equilibrium. A more stable (and costly) scheme is to do some streaming steps on the density field but not on the velocity field, before starting the actual simulations (see e.g. the Book of Guo and Shu).");
		enter_subsection("Iterative initialization stop condition");
		{
			declare_entry("Residual", "1e-7", dealii::Patterns::Double(1e-25),
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
		declare_entry("Number of time steps", "100000000",
				dealii::Patterns::Integer(1),
				"The maximum number of time steps.");
		declare_entry("Simulation end time", "100000000.0",
				dealii::Patterns::Double(0),
				"The end time of the simulation. "
						"Especially for adaptive time stepping schemes, number of steps is not an appropriate stop condition");
		declare_entry("Convergence threshold", "1e-30",
				dealii::Patterns::Double(),
				"The codes stops when the maximum velocity variation is below this threshold in 10 iterations.");
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
		declare_entry("Output checkpoint interval", "1000000000",
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
		enter_subsection("Turbulence Statistics");
		{
			declare_entry("Output global turbulence statistics?", "false",
					dealii::Patterns::Bool(),
					"Specifies if global turbulence statistics should be monitored.");
			declare_entry("Output turbulence statistics?", "false",
					dealii::Patterns::Bool(),
					"Specifies if turbulence statistics in slices should be monitored.");
			declare_entry("Wall normal direction", "1",
					dealii::Patterns::Integer(0, 3),
					"Convergence is monitored by putting out the turbulence statistics over planes that are parallel to the wall. The wall normal direction can be 0,1,2 for x,y,z, respectively.");
			declare_entry("Wall normal coordinates", "1e-1, 2e-1, 5e-1",
					dealii::Patterns::List(
							dealii::Patterns::Double(-1e10, 1e10)),
					"Convergence is monitored by putting out the turbulence statistics over planes that are parallel to the wall. This comma-separated list of decimal numbers specifies their wall-normal coordinates.");
		}
		leave_subsection();
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

	ParameterHandler::parse_input(filename);

} /* readFromTextFile */

/**
 * @short wrapper function for ParameterHandler::read_input_from_xml; directing cerr into a C++-Exception
 **/
void SolverConfiguration::readFromXMLFile(const std::string & filename) {

	// create stream
	std::ifstream xmlFile(filename);
	ParameterHandler::parse_input_from_xml(xmlFile);

} /* readFromXMLFile */

void SolverConfiguration::prepareOutputDirectory() {
	// make sure that the code is compiled with mpi
	assert(dealii::Utilities::MPI::job_supports_mpi());

	// Everything is checked by process 0.

	//  ((Using boost::filesystem provides a cross-platform solution))
	boost::filesystem::path outputDir(getOutputDirectory());

	// Make output directory
	if (is_MPI_rank_0()) {
		/// If not exists, try to create output directory
		boost::filesystem::path parentDir(outputDir.branch_path());
		if (not boost::filesystem::is_directory(parentDir)) {
			std::stringstream msg;
			msg << "You want to put your output directory into "
					<< parentDir.string()
					<< ", but this parent directory does not even exist.";
			throw ConfigurationException(msg.str());
		}
		try {
			//create_directory throws basic_filesystem_error<Path>, if fail (= no writing permissions)
			//returns false, if directory already existed
			boost::filesystem::create_directory(outputDir);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "You want to put your output directory into "
					<< parentDir.string()
					<< ", but you seem to have no writing permissions.";
			throw ConfigurationException(msg.str());
		}
		// Postcondition: directory exists
		try {
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
		if ((0 == getRestartAtIteration())
				and (not boost::filesystem::is_empty(outputDir))) {
			if (isUserInteraction()) {
				// Request user input
				pout
						<< "'Restart at checkpoint' is disabled, but Output directory is not empty. "
								"The simulation might overwrite old data. Do you really want to continue?"
								"If you are running your simulation in a parallel environment, you might want to "
								"switch user interaction off (which can be done by the corresponding option in SolverConfiguration"
								"or in the parameter file)." << endl;
				size_t yes1_or_no2 = 0; // = 1 for yes; = 2 for no
				string input = "";
				for (size_t i = 0; true; i++) {
					pout << "Please enter 'y' or 'n':" << endl;
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
					pout << "Your input was not understood. ";
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
			} /* if/else user interaction */
		} /* if not is restart at checkpoint */
		// Check writing permissions in directory (for every MPI process)
		try {
			/// try to create a single file
			std::ofstream filestream;
			std::stringstream filename;
			filename << "testfile_process."
					<< dealii::Utilities::MPI::this_mpi_process(
					MPI_COMM_WORLD) << ".txt";
			filestream.open(
					(outputDir / filename.str().c_str()).string().c_str());
			filestream << " ";
			filestream.close();
			boost::filesystem::remove(
					(outputDir / filename.str().c_str()).string().c_str());
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "Process " << dealii::Utilities::MPI::this_mpi_process(
			MPI_COMM_WORLD) << " running on host "
					<< dealii::Utilities::System::get_hostname()
					<< " does not have writing access to your Output directory "
					<< outputDir.string();
			throw ConfigurationException(msg.str());
		} /* catch */
		// make checkpoint dir
		try {
			boost::filesystem::path checkpoint_dir(outputDir / "checkpoint");
			boost::filesystem::create_directory(checkpoint_dir);
		} catch (std::exception& e) {
			std::stringstream msg;
			msg << "You want to put your checkpoint directory into "
					<< outputDir.string()
					<< ", but you seem to have no writing permissions.";
			throw ConfigurationException(msg.str());
		}
	} /* if is_MPI_rank_0() */
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
						"time integration schemes, you will have to set Time integrator to 'OTHER'."
				<< endl;
	}
	// check if a forcing scheme is set to perform pseudopotential collisions
	if (BGK_MULTIPHASE == getCollisionScheme()) {
		if (NO_FORCING == getForcingScheme()) {
			std::stringstream msg1;
			msg1
					<< "If you use the BGK Multiphase model you will need to specify a forcing scheme. Found 'No Forcing'.";
			LOG(ERROR) << msg1.str().c_str() << endl;
			throw ConfigurationException(msg1.str());
		}
		LOG(WARNING)
				<< "BGK_MULTIPHASE does not work properly. The finite element discretization "
						"does not provide gradients with sufficient isotropy"
				<< endl;
	}
	// not all forcing schemes work properly
	if (NO_FORCING != getForcingScheme()) {
		if (Stencil_D3Q15 == getStencil()) {
			LOG(ERROR)
					<< "NATriuM's Stencil_D3Q15 does not support forcing schemes, so far."
					<< endl;
		}
		if (Stencil_D3Q27 == getStencil()) {
			LOG(ERROR)
					<< "NATriuM's Stencil_D3Q27 does not support forcing schemes, so far."
					<< endl;
		}
		if (Stencil_D3Q19 == getStencil()) {
			if (SHIFTING_VELOCITY != getForcingScheme())
				LOG(ERROR)
						<< "NATriuM's Stencil_D3Q19 works only with the shifting velocity forcing scheme."
						<< endl;
		}
	}

	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1) {
		if (this->isUserInteraction()) {
			std::stringstream msg;
			msg
					<< "UserInteraction is switched on in the solver configuration. "
					<< "As that might be very dangerous when working on multiple MPI processes, "
					<< "I am throwing an exception (Safety first!) Please change your solver "
					<< "configuration (either in the parameter file or in natrium::SolverConfiguration.)"
					<< "It might be a good idea to make sure (manually) that you are not overwriting anything. ";
			LOG(ERROR) << msg.str() << endl;
			throw ConfigurationException(msg.str());
		}
	}
}

} /* namespace natrium */
