/**
 * @file step-2.cpp
 * @short Second tutorial:  Poiseuille Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

#include "TurbulentChannelFlow3D.h"
#include "FinalChannelStatistics.h"

#include "deal.II/numerics/data_out.h"
#include "deal.II/grid/grid_out.h"

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

	//pout << "Usage: ./turbulent-channel3D <is_restarted (false: 0, true: 1) <refinementLevel> <p> <filterID (no: 0, exp: 1, new: 2)>" << endl;
	pout << "Starting NATriuM step-turbulent-channel..." << endl;

	//**** User Input ****
	/**
	 * @ Simulation setup due to KMM (AIP 1998)
	 *  					| case #1			| case #2			| case #3
	 *  ==========================================================================
	 *  ReTau 				| 180				| 395				| 590
	 *  u_cl				| 10				| 10				| 10
	 *  uCl2uTauRatio		| 18.3				| 20.1332			| 21.2631
	 *
	 *  height				| 1					| 1					| 1
	 *  length				| 4					| 2					| 2
	 *  width				| 4/3				| 1					| 1
	 *
	 *  repetitions.at(0)	| 1					| 4					| 3
	 *  repetitions.at(1)	| 1					| 3					| 2
	 *  repetitions.at(2)	| 1					| 3					| 3
	 *
	 *  Ma					| 0.05;				| 0.05				| 0.05
	 */

	// Flow variables
	const double CFL = 2.0; 
	const double ReTau = atof(argv[1]);
	const double u_cl = atof(argv[2]);
	const double uCl2uTauRatio = atof(argv[3]);

	// Computational domain
	const double height = atof(argv[4]);
	const double length = atof(argv[5]) * M_PI * height / 2;
	const double width = atof(argv[6]) * M_PI * height / 2;

	// Grid resolution
	std::vector<unsigned int> repetitions(3);
	repetitions.at(0) = atoi(argv[7]);
	repetitions.at(1) = atoi(argv[8]);
	repetitions.at(2) = atoi(argv[9]);

	const double Ma = atof(argv[10]);// lower Ma => reduction in numerical compressibility

	const int refinementLevel = atoi(argv[11]);
	const int orderOfFiniteElement = atoi(argv[12]);
	const int filterID = atoi(argv[13]);

	int restart_iteration = atoi(argv[14]);
	bool is_periodic = true;

	// Turbulence statistics
	int noSamplePoints = 3;
	std::vector<double> samplePointCoordinates(noSamplePoints);

	if (ReTau == 180) {
		samplePointCoordinates[0] = 8. / 16;
		samplePointCoordinates[1] = 4. / 16;
		samplePointCoordinates[2] = 1. / 16;
	} else if (ReTau == 395) {
		samplePointCoordinates[0] = 12. / 24;
		samplePointCoordinates[1] = 6. / 24;
		samplePointCoordinates[2] = 1. / 24;
	}


	for (int i = 0; i < noSamplePoints; i++) {
		samplePointCoordinates[i] = 0.5 * height
				* (1 - cos( M_PI / height * samplePointCoordinates[i]));
	}

	//**** Calculated ****
	// approximate air viscosity at room temperature (275K): 1.3e-5 [m^2/s]
	double viscosity = u_cl * height / 2 / (ReTau * uCl2uTauRatio);
	//TODO: smooth increase of the inlet velocity until the initTime is reached
	//  	e.g. u_cl_init = u_cl*(F1B2 - F1B2 * cos(PI / (initTime * globalTimeStep)) ;

	const double scaling = sqrt(3) * u_cl / Ma;

	// Display user input
	pout << "============================================================="
			<< "\n" << " READ COMMAND LINE PARAMETERS " << "\n"
			<< "============================================================="
			<< "\n" << " |  Parameter \t\t\t\t| Value" << "\n"
			<< " +--------------------------------------+------------------- "
			<< "\n" << " |  Friction Reynolds number ReTau \t| " << ReTau
			<< "\n" << " |  Mean centerline velocity u_cl \t| " << u_cl << "\n"
			<< " |  Center line velocity to \t\t| " << "\n"
			<< " |  ... friction velocity ratio \t| " << uCl2uTauRatio << "\n"
			<< " |  \t\t\t\t\t| " << "\n" << " |  Channel height \t\t\t| "
			<< height << "\n" << " |  Channel length \t\t\t| " << length << "\n"
			<< " |  Channel width \t\t\t| " << width << "\n"
			<< " |  \t\t\t\t\t| " << "\n" << " |  Mach number Ma \t\t\t| " << Ma
			<< "\n" << " |  Repetitions at x \t\t\t| " << repetitions.at(0)
			<< "\n" << " |  Repetitions at y \t\t\t| " << repetitions.at(1)
			<< "\n" << " |  Repetitions at z \t\t\t| " << repetitions.at(2)
			<< "\n" << " |  Refinement level N \t\t\t| " << refinementLevel
			<< "\n" << " |  Order of finite element p \t\t| "
			<< orderOfFiniteElement << endl;

	///**** Create CFD problem ****
	boost::shared_ptr<TurbulentChannelFlow3D> channel3D = boost::make_shared<
			TurbulentChannelFlow3D>(viscosity, refinementLevel, repetitions,
			ReTau, u_cl, height, length, width, orderOfFiniteElement,
			is_periodic);
	//channel3D->refineAndTransform();

	//std::ofstream out_file("/tmp/grid_out.vtk");
	//dealii::GridOut().write_vtk(*channel3D->getMesh(), out_file);
	//out_file.close();
	//return 0;


	//viscosity = 0.5*dt*scaling*scaling/3.; //u_bulk * height / Re;
	//poiseuille2D->setViscosity(viscosity);
	//poiseuille2D->getExternalForce()->scale(viscosity);

	/// setup configuration
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/turbulent-channel3D/Re" << ReTau
			<< "-N" << refinementLevel << "-p" << orderOfFiniteElement
			<< "-filt" << filterID;

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	//configuration->setSwitchOutputOff(true);
	configuration->setOutputDirectory(dirName.str());
        cout << "Restart iteration " << restart_iteration << endl;
	configuration->setRestartAtIteration(restart_iteration);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(100);
	configuration->setOutputCheckpointInterval(40000);
	configuration->setOutputSolutionInterval(1000);
	configuration->setCommandLineVerbosity(WELCOME);
	configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setCFL(CFL);
	configuration->setForcingScheme(SHIFTING_VELOCITY);
	configuration->setStencil(Stencil_D3Q19);


	if (filterID == 1) {
		configuration->setFiltering(true);
		configuration->setFilteringScheme(EXPONENTIAL_FILTER);
	} else if (filterID == 2) {
		configuration->setFiltering(true);
		configuration->setFilteringScheme(NEW_FILTER);
	}

	configuration->setOutputTurbulenceStatistics(true);
	configuration->setWallNormalDirection(1);
	configuration->setWallNormalCoordinates(samplePointCoordinates);

	configuration->setTimeIntegrator(EXPONENTIAL);
	//configuration->setDealIntegrator(SDIRK_TWO_STAGES);
	//configuration->setInitializationScheme(ITERATIVE);
	//configuration->setIterativeInitializationNumberOfIterations(100);
	//configuration->setIterativeInitializationResidual(1e-15);
	configuration->setConvergenceThreshold(1e-10);
	//configuration->setNumberOfTimeSteps(1);

	//configuration->setSimulationEndTime(); // unit [s]

	// ----------------------------------------------------------
	// create a separate object for the initial velocity function
	TurbulentChannelFlow3D::IncompressibleU test_velocity(channel3D.get());

	if (restart_iteration == 0) {
		// Divergence check
		pout << "**** Divergence check ****" << endl;
		//srand(1);
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
			div += ((f_h - f) / h);
			x_plus_h(0) = x(0);

			// dv / dy
			x_plus_h(1) = x(1) + h;
			f = test_velocity.value(x, 1);	// component 1 -> v
			f_h = test_velocity.value(x_plus_h, 1);
			div += ((f_h - f) / h);
			x_plus_h(1) = x(1);

			// dw / dz
			x_plus_h(2) = x(2) + h;
			f = test_velocity.value(x, 2);	// component 2 -> w
			f_h = test_velocity.value(x_plus_h, 2);
			div += ((f_h - f) / h);

			// check div small (could also be done with asserts)
			pout << "... div u at point " << i << ", y-coord " << x(2) << ": "
					<< div << endl;
		}
	}

	// make solver object and run simulation
	CFDSolver<3> solver(configuration, channel3D);

	if (restart_iteration == 0) {
		double utrp_max = channel3D.get()->getMaxUtrp();
		double utrp_inc_max = channel3D.get()->getMaxIncUtrp();
		double scalingFactor = utrp_max / utrp_inc_max;

		pout << " >>>> Max. Velocity Perturbation: " << utrp_max << endl;
		pout << " >>>> Max. Incompressible Velocity Perturbation:  "
				<< utrp_inc_max << endl;
		pout << " >>>> Scaling Factor:  " << scalingFactor << endl;

		solver.scaleVelocity(scalingFactor);
		solver.addToVelocity(
				boost::make_shared<TurbulentChannelFlow3D::MeanVelocityProfile>(
						channel3D.get()));
	} else {
		solver.appendDataProcessor(
				boost::make_shared<FinalChannelStatistics>(solver,
						configuration->getOutputDirectory()));
	}
	solver.run();

	pout << "Max Velocity  " << solver.getMaxVelocityNorm() << endl; //"   (laminar: "<< 1.5*u_bulk << ")" <<  endl;

	pout << "NATriuM step-turbulent-channel terminated." << endl;

	return 0;
}
