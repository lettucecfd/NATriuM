/**
 * @file step-2.cpp
 * @short Second tutorial:  Poiseuille Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>

#include "TurbulentChannelFlow3D.h"
#include "FinalChannelStatistics.h"
#include "AdaptiveForcing.h"

#include "deal.II/numerics/data_out.h"
#include "deal.II/grid/grid_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/D3Q19.h"

#include "natrium/problemdescription/ProblemDescription.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/ConfigNames.h"
#include "natrium/utilities/CommandLineParser.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

	MPIGuard::getInstance(argc, argv);

	// ========================================================================
	// READ COMMAND LINE PARAMETERS
	// ========================================================================

	CommandLineParser parser(argc, argv);
	parser.addDocumentationString("turbulent-channel",
			"Turbulent channel flow");
	parser.setPositionalArgument<int>("ref-level",
			"Global refinement level. Number of total grid points in direction i: 'ref-i' * (2 ^ 'ref-level')");
	parser.setArgument<int>("restart", "Restart at iteration ...", 0);
    parser.setArgument<int>("compressible", "Is a thermal/compressible LB model to be used?",0);
	parser.setArgument<double>("Re_tau",
			"Wall Reynolds number u_tau * delta / nu", 180.0);
	parser.setArgument<double>("U_cl", "Centerline velocity for initialization",
			10.0);
	parser.setArgument<double>("delta", "Channel half width", 1);
	parser.setArgument<double>("lx",
			"Streamwise length of computational domain normalized with delta",
			4 * M_PI);
	parser.setArgument<double>("lz",
			"Spanwise length of computational domain normalized with delta",
			2 * M_PI);
	parser.setArgument<double>("Re_cl",
			"Centerline Reynolds number used for initialization Re_c = U_cl * delta / nu",
			3300);
	parser.setArgument<double>("Ma", "Mach number U_cl/cs", 0.1);
    parser.setArgument<double>("grid-density", "The higher the value, the denser the grid at the walls", 0.8);
	parser.setArgument<int>("rep-x",
			"Number of repetitions in x-direction (to refine the grid in steps that are not 2^N).",
			1);
	parser.setArgument<int>("rep-y", "cf. rep-x", 1);
	parser.setArgument<int>("rep-z", "cf. rep-x", 1);
	parser.setFlag("test", "no simulation, just show message");
	try {
		parser.importOptions();
	} catch (HelpMessageStop&) {
		return 0;
	}

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

	const int restart = parser.getArgument<int>("restart");
	LOG(WELCOME) << "==================================================="
			<< endl << "=== Starting NATriuM step-turbulent-channel... ===="
			<< endl << "=== Restart iteration: " << restart << endl
			<< "===================================================" << endl;

	// ========================================================================
	// CHANNEL SETUP
	// ========================================================================
    const double gamma = 1.4;
	// Reynolds number
	const double Re_tau = parser.getArgument<double>("Re_tau");
	assert(Re_tau > 0);

	// Computational domain
	const double delta = parser.getArgument<double>("delta");
	const double height = 2.0 * delta;
	const double length = parser.getArgument<double>("lx") * delta;
	const double width = parser.getArgument<double>("lz") * delta;

	// Grid resolution
	const int ref_level = parser.getArgument<int>("ref-level");
	std::vector<unsigned int> repetitions(3);
	repetitions.at(0) = parser.getArgument<int>("rep-x");
	repetitions.at(1) = parser.getArgument<int>("rep-y");
	repetitions.at(2) = parser.getArgument<int>("rep-z");

	// calculate viscosity and scaling
	const double u_cl = parser.getArgument<double>("U_cl");
	const double Ma = parser.getArgument<double>("Ma");
	const double gridDensity = parser.getArgument<double>("grid-density");
	const double scaling = sqrt(3) * u_cl / (Ma*sqrt(gamma));

	const double viscosity = u_cl * delta / parser.getArgument<double>("Re_cl");
	const double utau = Re_tau * viscosity / delta;

	// make channel flow object
	boost::shared_ptr<TurbulentChannelFlow3D> channel3D = boost::make_shared<
			TurbulentChannelFlow3D>(viscosity,
                                    (size_t) parser.getArgument<int>("ref-level"), repetitions, Re_tau,
                                    u_cl, height, length, width, gridDensity);

	// Turbulence statistics
	// y-coordinates for output of RMS values in table (turbulence monitor)
	int n_rms_coords = 3;
	std::vector<double> rms_coords(n_rms_coords);
	rms_coords[0] = height / 2.0;
	rms_coords[1] = height / 8.0;
	rms_coords[2] = height / 32.0;
	TurbulentChannelFlow3D::UnstructuredGridFunc trafo(length, height, width);
	for (int i = 0; i < n_rms_coords; i++) {
		rms_coords.at(i) = trafo.trans(rms_coords.at(i));
	}
	double ymin = trafo.trans(height / repetitions.at(1) / pow(2, ref_level));
	double yplus = ymin / (viscosity / utau);

	// ========================================================================
	// SOLVER CONFIGURATION
	// ========================================================================

	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();
	configuration->setRestartAtIteration(restart);
	configuration->setUserInteraction(false);
	configuration->setOutputTableInterval(100);
	configuration->setOutputCheckpointInterval(5000);
	configuration->setOutputSolutionInterval(5000);
	configuration->setCommandLineVerbosity(WELCOME);
	configuration->setStencilScaling(scaling);
	configuration->setCommandLineVerbosity(ALL);
	configuration->setForcingScheme(EXACT_DIFFERENCE);
    configuration->setStencil(Stencil_D3Q45);
    configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
    configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
    //configuration->setPrandtlNumber(0.7);
    configuration->setSutherlandLaw();
	configuration->setOutputTurbulenceStatistics(true);
	configuration->setWallNormalDirection(1);
	configuration->setWallNormalCoordinates(rms_coords);
	parser.applyToSolverConfiguration(*configuration);
	std::stringstream dir;
	dir << getenv("NATRIUM_HOME") << "/KMM/Re" << Re_tau << "-ref" << ref_level
			<< "-p" << configuration->getSedgOrderOfFiniteElement() << "-Ma"
			<< Ma << "-cfl" << configuration->getCFL() << "-rep"
			<< repetitions.at(0) << repetitions.at(1) << repetitions.at(2) << "-sten"
            << configuration->getStencil() << "-gridDensity" << gridDensity;
	configuration->setOutputDirectory(dir.str());

	// ========================================================================
	// COMMAND LINE OUTPUT
	// ========================================================================
	LOG(WELCOME) << "          -----         " << endl
			<< "          -----         " << endl << "FLOW SETUP: " << endl
			<< "===================================================" << endl
			<< "Re_cl = u_cl * delta / nu   = " << u_cl << " * " << delta
			<< " / " << viscosity << " = " << u_cl * delta / viscosity << endl
			<< "u_tau = Re_tau * nu / delta = " << Re_tau << " * " << viscosity
			<< " / " << delta << " = " << Re_tau * viscosity / delta << endl
			<< "F     = rho * utau^2 / delta = 1.0 * " << utau << "^2" << " / "
			<< delta << " = " << utau * utau / delta << endl
			<< "          -----         " << endl << "          -----        "
			<< endl;
	const double dxplus = length / repetitions.at(0) / pow(2, ref_level)
			/ (viscosity / utau);
	const double dzplus = width / repetitions.at(2) / pow(2, ref_level)
			/ (viscosity / utau);
	const double p = configuration->getSedgOrderOfFiniteElement();
	LOG(WELCOME) << "CHANNEL SETUP: " << endl
			<< "===================================================" << endl
			<< "Dimensions:    " << length << " x " << height << " x " << width
			<< endl << "Grid:          " << repetitions.at(0) << " x "
			<< repetitions.at(1) << " x " << repetitions.at(2)
			<< " blocks with 8^" << ref_level << " cells each " << endl
			<< "#Cells:        " << int(repetitions.at(0) * pow(2, ref_level))
			<< " x " << int(repetitions.at(1) * pow(2, ref_level)) << " x "
			<< int(repetitions.at(2) * pow(2, ref_level)) << " = "
			<< int(
					repetitions.at(0) * repetitions.at(1) * repetitions.at(2)
							* pow(2, 3 * ref_level)) << endl << "#Points:       "
			<< int(repetitions.at(0) * pow(2, ref_level) * p) << " x "
			<< int(repetitions.at(1) * pow(2, ref_level) * p) << " x "
			<< int(repetitions.at(2) * pow(2, ref_level) * p) << " = "
			<< int(
					repetitions.at(0) * repetitions.at(1) * repetitions.at(2)
							* pow(2, 3 * ref_level) * p * p * p) << endl
			<< "y+ (wrt. cells): "  << yplus << "   dx+ = " << dxplus << ", "
			<< "dz+ = " << dzplus << endl << "          -----         " << endl
			<< "          -----         " << endl;
	const double dt = configuration->getCFL() / (p * p)
			/ (sqrt(2) * scaling) * ymin;
	LOG(WELCOME) << "==================================================="
			<< endl << "               dt  = " << dt << endl << "               dt+ = "
			<< dt/(viscosity/utau/utau) << endl << "    u_tau cross time = "
			<< length / utau << " = "
			<< int(length / utau / dt) << " steps" << endl
			<< "===================================================" << endl
			<< endl;

	if (parser.hasArgument("test")) {
		return 0;
	}

	// ----------------------------------------------------------
	// create a separate object for the initial velocity function
	// TurbulentChannelFlow3D::IncompressibleU test_velocity(channel3D.get());

	/*if (restart == 0) {
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
	 } */

	// ========================================================================
	// CREATE SOLVER AND RUN SIMULATION
	// ========================================================================
	// make solver object and run simulation
	CompressibleCFDSolver<3> solver(configuration, channel3D);

	if (restart == 0) {
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
	} /* else {
		// append data processor only in case of restart
		solver.appendDataProcessor(
				boost::make_shared<FinalChannelStatistics>(solver,
						configuration->getOutputDirectory()));
	}*/
    solver.appendDataProcessor(
            boost::make_shared<AdaptiveForcing>(solver, configuration->getOutputDirectory(), 1.0)); 
	solver.run();

	pout << "Max Velocity  " << solver.getMaxVelocityNorm() << endl;

	pout << "NATriuM step-turbulent-channel terminated." << endl;

	return 0;
}
