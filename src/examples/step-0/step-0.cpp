/**
 * @file step-0.cpp
 * @example step-0.cpp
 * @short
 * @section sec_intro Basis tutorial: Lid-driven cavity in 2D
 *    This example is thoroughly commented to introduce you to the key concepts in
 *    fluid flow simulations with NATriuM. There are three main classes that are
 *    used in a NATriuM simulation:
 *    - SolverConfiguration: Contains all the simulation preferences.
 *    - ProblemDescription: Contains the computation grid, initial conditions and
 *      boundary conditions. ProblemDescription is an abstract class with the purely
 *      virtual functions applyInitialDensities and applyInitialVelocities.
 *      When you use NATriuM as a library, this class has to be overriden by a problem specific
 *      child class (here the class LidDrivenCavity2D), which implements the virtual
 *      functions.
 *    - CFDSolver: The CFDSolver is the central class in the NATriuM framework.
 *      It creates all necessary classes and runs the simulation.
 *
 *    We have the following includes that are part of every simulation script in NATriuM:
 *    First, we need stdlib.h to access the environment variable NATRIUM_HOME
 *    and ssteam to concatenate strings to build the path of the output directory.
 *    Second, the three classes mentioned above are implemented in CFDSolver.h,
 *    SolverConfiguration.h, and ProblemDescription.h.
 *    Third, BasicNames.h contains some general typedefs used throughout NATriuM and
 *    automatically enables using certain std and boost types and functions (e.g. shared_ptr).
 *    It has to be included in basically every file that is part of or uses NATriuM.
 *    Finally, we include LidDrivenCavity2D.h. It contains the class LidDrivenCavity2D,
 *    which defines the concrete flow we want to simulate.
 *          @snippet step-0.cpp Includes
 *
 *    We include the namespace natrium to access all classes and functions directly.
 *    		@snippet step-0.cpp Namespace
 *
 *    The Main function of the program begins with defining the Reynolds and Mach number.
 *    Note that in the LBM context the Mach number must not exceed a limit of approximately
 *    0.3. On the other hand the time step is proportional to the Mach number, meaning
 *    that for Ma << 0.1 the runtime of the code will grow. Usually, Ma is set roughly to 0.1.
 *    (On the contrary, the Reynolds number does not influence the time step size directly.)
 *    		@snippet step-0.cpp Main function
 *
 *    Next, the discretization is defined. As the grid generation in the easier geometries
 *    is done inside deal.II, we need to specify the refinement level. Refinement level 4 will
 *    mean that our quadratic domain will have 2^4 * 2^4 =  256 cells. See LidDrivenCavity2D::makeGrid()
 *    for details on the grid generation. The next thing to specify is the order of the
 *    finite element, often referred to as p. p=1 stands for linear shape functions (the
 *    finite elments), while greater values of p stand for higher order polynomials.
 *    The choice of p has a tremendous impact on the accuracy and runtime of the program.
 *    Great values of p can lead to very accurate solutions (sometimes even exponential
 *    convergence in space!) but they require tiny timesteps (\f[ \delta_t \sim  \frac{1}{(p+1)^2} \f]).
 *    To chose p=5 is often a reasonable compromise.
 *    Finally, we have to chose the time step. Small timesteps lead to long simulations, while big
 *    timesteps render the simulation unstable. When we use explicit time integration, the suitable time
 *    step depends only on the CFL condition. There are functions implemented in NATriuM to automate
 *    theses calculations. This is left for later tutorials and we simply set dt=0.0001.
 *    		@snippet step-0.cpp Discretization
 *
 *	  Now, we have to specify the characteristic velocity and kinematic viscosity of the flow. The setting
 *	  of the characteristic velocity depends on the fact, that we fix the speed of sound to
 *	  \f[ c_s = \frac{1}{\sqrt{3}} \f] and normalize the time scale with respect to that speed.
 *	  This is a standard approach in Lattice Boltzmann. However, NATriuM also provides the possibility
 *	  to use physical units, which has to do with the scaling of the difference stencil is covered in later tutorials.
 *	  The viscosity is determined by the Reynolds number, characteristic velocity and characteristic length L = 1.0.
 *    		@snippet step-0.cpp Definition
 *
 *    The next step is to create the configuration object and feed it with the information we have just discussed.
 *    There are a lot more things to set, but we take the default settings for the rest. The the output frequency
 *    and number of time steps are set manually. Additionally, we disable the restart option, which forces the
 *    simulation to start at t=0.
 *    		@snippet step-0.cpp Configuration
 *
 *	  Similarly, we create the flow object. It is inherited of ProblemDescription and needed to be assigned to
 *	  a shared_ptr<ProblemDescription>.
 *    		@snippet step-0.cpp Problem
 *
 *	  Finally, we are ready to create the CFD solver and run the simulation.
 *    		@snippet step-0.cpp Solver
 *
 *
 *
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */


//! [Includes]
#include <stdlib.h>
#include <sstream>

#include "solver/CFDSolver.h"
#include "solver/SolverConfiguration.h"

#include "problemdescription/ProblemDescription.h"

#include "utilities/BasicNames.h"

#include "LidDrivenCavity2D.h"
//! [Includes]

//! [Namespace]
using namespace natrium;
//! [Namespace]



//! [Main function]
int main() {

	cout << "Starting NATriuM step-0..." << endl;

	// set Reynolds and Mach number
	const double Re = 1000;
	const double Ma = 0.1;
	//! [Main function]

	//! [Discretization]
	// set spatial discretization
	size_t refinementLevel = 4;
	size_t orderOfFiniteElement = 5;

	// set temporal discretization
	double dt = 0.0001;
	//! [Discretization]

	//! [Definition]
	// set Problem so that the right Re and Ma are achieved
	const double U = 1 / sqrt(3) * Ma;
	const double viscosity = U / Re; // (because L = 1)

	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/step-0";
	//! [Definition]

	//! [Configuration]
	shared_ptr<SolverConfiguration> configuration = make_shared<
			SolverConfiguration>();
	configuration->setOutputDirectory(dirname.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setOutputCheckpointInterval(1000);
	configuration->setOutputSolutionInterval(100);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setTimeStepSize(dt);
	configuration->setNumberOfTimeSteps(200000);
	//! [Configuration]

	//! [Problem]
	// make problem and solver objects
	shared_ptr<LidDrivenCavity2D> lidDrivenCavity = make_shared<
			LidDrivenCavity2D>(U, viscosity, refinementLevel);
	shared_ptr<ProblemDescription<2> > ldCavityProblem = lidDrivenCavity;
	//! [Problem]

	//! [Solver]
	CFDSolver<2> solver(configuration, ldCavityProblem);
	solver.run();

	cout << "NATriuM step-0 terminated." << endl;

	return 0;
}

//! [Solver]
