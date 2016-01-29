/**
 * @file step-2.cpp
 * @short Second tutorial:  Poiseuille Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

#include "TurbulentChannelFlow3D.h"

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/D3Q19.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"


using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	pout << "Starting NATriuM step-turbulent-channel..." << endl;

	//**** User Input ****
	const double CFL 					= 0.4;
	const double ReTau 					= 180;
	const double ReCl 					= 3300;
	const double Re_bulk 				= 5600;
	const double u_bulk					= 10;
	const double height 				= 1;
	const double length 				= 2 * 3 * height; 	// PI = 3;
	const double width 					= 3 * height;		// PI = 3;
	const double orderOfFiniteElement 	= 2;
	const double Ma 					= atof(argv[1]);
	const double refinement_level 		= atoi(argv[2]);
	bool is_periodic 					= true;

	//**** Calculated ****
	// approximate air viscosity at room temperature (275K): 1.3e-5 [m^2/s]
	double viscosity  = u_bulk * height / Re_bulk;

	//TODO: smooth increase of the inlet velocity until the initTime is reached
	//  	e.g. u_cl_init = u_cl*(F1B2 - F1B2 * cos(PI / (initTime * globalTimeStep)) ;

	///**** Create CFD problem ****
	const double scaling = sqrt(3) * 1.16 * u_bulk / Ma; //TODO: not laminar!
	boost::shared_ptr<TurbulentChannelFlow3D> channel3D =
			boost::make_shared<TurbulentChannelFlow3D>(viscosity, refinement_level, ReTau, ReCl,
					u_bulk, height, length, width, orderOfFiniteElement, is_periodic);
	const double dt = CFDSolverUtilities::calculateTimestep<3>(
			*channel3D->getMesh(), orderOfFiniteElement, D3Q19(scaling), CFL);

	//viscosity = 0.5*dt*scaling*scaling/3.; //u_bulk * height / Re;
	//poiseuille2D->setViscosity(viscosity);
	//poiseuille2D->getExternalForce()->scale(viscosity);


	// -----------------------------------------------------------------------------------------------------------------------
	// create a separate object for the initial velocity function (Constructor has to get the "flow" object or a pointer to it or something)
	TurbulentChannelFlow3D::IncompressibleU test_velocity(channel3D.get());

	// Divergence check
	cout << "**** Divergence check ****" << endl;
	srand(1);
	for (size_t i = 0; i < 30; i++) {

		double div = 0.0;

		// create random point in the flow domain and calculate f
		dealii::Point<3> x;
		x(0) = (double) random() / RAND_MAX * length;
		x(1) = (double) random() / RAND_MAX * width;
		x(2) = (double) random() / RAND_MAX * height;

		// Calculate div(U):

		// increment
		double h = 1e-6;
		dealii::Point<3> x_plus_h(x);

		// du / dx
		x_plus_h(0) = x(0) + h;
		double f = test_velocity.value(x, 0); // component 0 -> u
		double f_h = test_velocity.value(x_plus_h, 0);
		div += ( (f_h - f) / h );
		x_plus_h(0) = x(0);

		// dv / dy
		x_plus_h(1) = x(1) + h;
		f = test_velocity.value(x, 1);	// component 1 -> v
		f_h = test_velocity.value(x_plus_h, 1);
		div += ( (f_h - f) / h );
		x_plus_h(1) = x(1);


		// dw / dz
		x_plus_h(2) = x(2) + h;
		f = test_velocity.value(x, 2);	// component 2 -> w
		f_h = test_velocity.value(x_plus_h, 2);
		div += ( (f_h - f) / h );

		// check div small (could also be done with asserts)
		pout << "... div u at point " << i << ", y-coord " << x(2) << ": "<< div << endl;
	}

	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/turbulent-channel3D";
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
	configuration->setRestartAtLastCheckpoint(false);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(100);
	configuration->setOutputCheckpointInterval(100000000);
	configuration->setOutputSolutionInterval(100);
	configuration->setCommandLineVerbosity(WELCOME);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setTimeStepSize(dt);
	configuration->setStencil(Stencil_D3Q19);
	//configuration->setTimeIntegrator(OTHER);
	//configuration->setDealIntegrator(CRANK_NICOLSON);

	//configuration->setInitializationScheme(ITERATIVE);
	//configuration->setIterativeInitializationNumberOfIterations(100);
	//configuration->setIterativeInitializationResidual(1e-15);

	configuration->setConvergenceThreshold(1e-10);
	//configuration->setNumberOfTimeSteps(100);
	//configuration->setSimulationEndTime(); // unit [s]

	// make solver object and run simulation
	CFDSolver<3> solver(configuration, channel3D);

	// DEBUG
//	TurbulentChannelFlow3D::InitialVelocity initUtrp(channel3D.get());
//	TurbulentChannelFlow3D::IncompressibleU initIncUtrp(channel3D.get());
//
//	double utrp_max = initUtrp.getMaxUtrp();
//	double utrp_inc_max = initIncUtrp.getMaxIncUtrp();

	//TurbulentChannelFlow3D maxUtrp(channel3D.get());

	double utrp_max = channel3D.get()->getMaxUtrp();
	double utrp_inc_max = channel3D.get()->getMaxIncUtrp();
	double scalingFactor = utrp_max / utrp_inc_max;

	pout << " >>>> Max. Velocity Perturbation: " << utrp_max << endl;
	pout << " >>>> Max. Incompressible Velocity Perturbation:  " << utrp_inc_max << endl;
	pout << " >>>> Scaling Factor:  " << scalingFactor << endl;

	solver.scaleVelocity(scalingFactor);
	solver.addToVelocity(boost::make_shared<TurbulentChannelFlow3D::MeanVelocityProfile>( channel3D.get() ) );
	solver.run();

	pout << "Max Velocity  " <<
			solver.getMaxVelocityNorm() << "   (laminar: "<< 1.5*u_bulk << ")" <<  endl;

	pout << "NATriuM step-turbulent-channel terminated." << endl;

	return 0;
}
