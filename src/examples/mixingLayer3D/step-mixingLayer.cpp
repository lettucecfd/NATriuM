//
// Created by dwilde3m on 02.12.21.
//
/*
 * step-mixingLayer.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */
#include <fstream>
//#include <time.h>
#include <stdlib.h>
//#include "deal.II/numerics/data_out.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/stencils/Stencil.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/CommandLineParser.h"
#include "MixingLayer3D.h"
#include "ShearLayerStats.h"
#include <unistd.h>
#include <mpi.h>


using namespace natrium;

double shearLayerThickness = 0.093;
//int dft_points = 20;

// Main function
int main(int argc, char** argv) {
    MPIGuard::getInstance(argc, argv);
    pout << "Starting NATriuM step-mixingLayer ..." << endl;
    /////////////////////////////////////////////////
    // read from command line
    //////////////////////////////////////////////////
    CommandLineParser parser(argc, argv);
    parser.setArgument<int>("Re", "Reynolds number 1/nu", 800);
    // TODO: Convective Mach number Mc=(U1-Uc)/c1, Uc=(U1c2+U2c1)/(c1+c2)
    //  Uc is the convective velocity of the large structures, U\
    //  and U2 are the freestream velocities, and c{ and c2 are the
    //  freestream sound speeds.
    // Set to 0.3, 0.7, 0.9, 1.0, 1.2
    parser.setArgument<double>("Ma", "Mach number", 0.3);
    parser.setArgument<double>("time", "simulation time (s)", 15);
    parser.setArgument<double>("randuscaling", "factor to scale random velocity field", 1);
    parser.setArgument<double>("uscaling", "factor to scale U1, i.e. deltaUx", 1);
    parser.setArgument<double>("CFL", "CFL number", 0.4);
    parser.setArgument<int>("nout", "output vtk every nout steps", 1000);
    parser.setArgument<int>("nstats", "output stats every nstats steps", 20);
    parser.setArgument<int>("squash", "squash grid towards centre", 0);
    parser.setArgument<int>("print", "print calculations of initial velocity", 0);
    parser.setArgument<int>("recalculate", "recalculate initial velocity", 0);
    parser.setArgument<string>("meshname", "name of the mesh file (shearlayer_*.txt)", "final_small");
    parser.setArgument<string>("randuname", "name of the initial velocity file (random_u_*.txt)", "k048_half");
    parser.setArgument<int>("order", "order of finite elements", 3);
    parser.setArgument<int>("ref-level", "Refinement level of the computation grid.", 0);
    parser.setArgument<int>("grid-repetitions",
                            "Number of grid cells along each axis before global refinement; "
                            "to produce grids with refinements that are not powers of two.", 3);
    parser.setArgument<int>("restart", "Restart at iteration ...", 0);

    try { parser.importOptions();
    } catch (HelpMessageStop&) { return 0;
    }
    auto meshname = parser.getArgument<string>("meshname");
    auto randuname = parser.getArgument<string>("randuname");
    double randuscaling = parser.getArgument<double>("randuscaling");
    double uscaling = parser.getArgument<double>("uscaling");
    double Re = parser.getArgument<int>("Re");
    double refinement_level = parser.getArgument<int>("ref-level");
    long nout = parser.getArgument<int>("nout");
    auto time = parser.getArgument<double>("time");
    const int restart = parser.getArgument<int>("restart");
    if (restart > 0) {
        LOG(WELCOME) << "==================================================="
                     << endl << "=== Starting NATriuM step-mixingLayer... ===="
                     << endl << "=== Restart iteration: " << restart << endl
                     << "===================================================" << endl;
    }

    /////////////////////////////////////////////////
    // set parameters, set up configuration object
    //////////////////////////////////////////////////
    // im Paper von Gassner und Beck ist U = 1/2pi definiert !!!!!
    // Aber ihre zeitangaben beziehen sich auf U = 1, wie bei Brachet (1991)
    // hier simulieren wir jetzt  U = 1 (ist im TGV3D modul sowieso nur so definiert)
    // und mit der Reynoldsähnlichkeit ist das kein Problem, da wir gegenüber Gassner
    // und Beck ja auch die Viskosität verändern
    // Nach einem Blick in van Rees et.al. (2011) und einer erfolgreichen Simulation in Palabos:
    // (Hier war tau = 3*nu_LB + 0.5, nu_LB = U_lattice * (N/2pi) / Re, dt = 2pi/N* U_lattice
    // Re = 1/nu,L=2pi, U = 1 und D = [0,2pi*L]^3
    const double U = 1;
    //const double L = 2 * M_PI;
    const double viscosity = 1.0 / Re;
    const double Ma = parser.getArgument<double>("Ma")*sqrt(1.4);
    const auto cfl = parser.getArgument<double>("CFL");
    const double cs = U / Ma;

    // chose scaling so that the right Ma-number is achieved
    const double reference_temperature = 1; //0.85;
    const double gamma = 1;//1.4;
    const double scaling = sqrt(3) * U / (Ma*sqrt(gamma*reference_temperature));
//    const double scaling = sqrt(3) * cs; // choose different? -> stencil larger/smaller -> from turb. channel

    // setup configuration
    boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<SolverConfiguration>();
    if (restart > 0) configuration->setRestartAtIteration(restart);
    configuration->setUserInteraction(false);
    configuration->setOutputCheckpointInterval(nout*100);
    configuration->setOutputSolutionInterval(nout);
    configuration->setSimulationEndTime(time);
    configuration->setOutputGlobalTurbulenceStatistics(true);
    configuration->setOutputCompressibleTurbulenceStatistics(true);
    configuration->setOutputShearLayerStatistics(true);
    configuration->setOutputShearLayerInterval(parser.getArgument<int>("nstats"));
    configuration->setStencilScaling(scaling);
    configuration->setStencil(Stencil_D3Q45);
    configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
    configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
    configuration->setHeatCapacityRatioGamma(1.4);
    configuration->setPrandtlNumber(0.71);
    configuration->setSedgOrderOfFiniteElement(parser.getArgument<int>("order")); // TODO: set to 4
    configuration->setCFL(cfl); // TODO: should be 0.4<CFL<2
//    configuration->setInitializationScheme(COMPRESSIBLE_ITERATIVE);

    parser.applyToSolverConfiguration(*configuration);

    // standard output dir
    string m_dirname;
    if (not parser.hasArgument("output-dir")){
        std::stringstream dirName;
        dirName << getenv("NATRIUM_HOME");
        dirName << "/step-mixingLayer/Re" << Re
                << "-Ma" << floor(Ma*1000)/1000
                << "-ref" << refinement_level
                << "-p" << configuration->getSedgOrderOfFiniteElement()
                << "-mesh" << meshname
                << "-randu" << randuname << "x" << floor(randuscaling*1000)/1000
                << "-uscale" << uscaling;
//        dirName << "-coll" << static_cast<int>(configuration->getCollisionScheme())
//                << "-sl" << static_cast<int>(configuration->getAdvectionScheme())
        if (configuration->getAdvectionScheme() != SEMI_LAGRANGIAN)
            dirName << "-int" << static_cast<int>(configuration->getTimeIntegrator()) << "_" << static_cast<int>(configuration->getDealIntegrator());
        dirName << "-CFL" << configuration->getCFL();
//        dirName << "-sten" << static_cast<int>(configuration->getStencil());
//        if (configuration->isFiltering()) (dirName << "-filt" << static_cast<int>(configuration->getFilteringScheme()) << "by_max_degree");
        if (configuration->getRegularizationScheme() != NO_REGULARIZATION)
            dirName << "-reg" << static_cast<int>(configuration->getRegularizationScheme());
//        if (configuration->getEquilibriumScheme()!= BGK_EQUILIBRIUM) (dirName << "-equili" << static_cast<int>(configuration->getEquilibriumScheme()));
        if (configuration->getCollisionScheme() == MRT_STANDARD) {
            dirName << "-mrt" << static_cast<int>(configuration->getMRTBasis());
        }
        if (configuration->getCollisionScheme() == MRT_STANDARD) {
            dirName << "-relax" << static_cast<int>(configuration->getMRTRelaxationTimes());
        }
        configuration->setOutputDirectory(dirName.str());
        m_dirname = dirName.str();
    } else {
        m_dirname = configuration->getOutputDirectory();
    }

    // ========================================================================
    // COMMAND LINE OUTPUT
    // ========================================================================
//    const double p = configuration->getSedgOrderOfFiniteElement();
//    const double dt = configuration->getCFL() / (p * p) / (sqrt(2) * scaling) * ymin;
//    const double dxplus = length / repetitions.at(0) / pow(2, ref_level)
//                          / (viscosity / utau);
//    const double dzplus = width / repetitions.at(2) / pow(2, ref_level)
//                          / (viscosity / utau);
//    LOG(WELCOME) << "          -----         " << endl
//                    << "          -----         " << endl << "FLOW SETUP: " << endl
//                    << "===================================================" << endl
//                    << "Re_cl = u_cl * delta / nu   = " << u_cl << " * " << delta
//                    << " / " << viscosity << " = " << u_cl * delta / viscosity << endl
//                    << "u_tau = Re_tau * nu / delta = " << Re_tau << " * " << viscosity
//                    << " / " << delta << " = " << Re_tau * viscosity / delta << endl
//                    << "F     = rho * utau^2 / delta = 1.0 * " << utau << "^2" << " / "
//                    << delta << " = " << utau * utau / delta << endl
//                    << "          -----         " << endl << "          -----        "
//                    << endl
//                    << "CHANNEL SETUP: " << endl
//                    << "===================================================" << endl
//                    << "Dimensions:    " << lx << " x " << height << " x " << width
//                    << endl << "Grid:          " << repetitions.at(0) << " x "
//                    << repetitions.at(1) << " x " << repetitions.at(2)
//                    << " blocks with 8^" << ref_level << " cells each " << endl
//                    << "#Cells:        " << int(repetitions.at(0) * pow(2, ref_level))
//                    << " x " << int(repetitions.at(1) * pow(2, ref_level)) << " x "
//                    << int(repetitions.at(2) * pow(2, ref_level)) << " = "
//                    << int(
//                         repetitions.at(0) * repetitions.at(1) * repetitions.at(2)
//                         * pow(2, 3 * ref_level)) << endl << "#Points:       "
//                    << int(repetitions.at(0) * pow(2, ref_level) * p) << " x "
//                    << int(repetitions.at(1) * pow(2, ref_level) * p) << " x "
//                    << int(repetitions.at(2) * pow(2, ref_level) * p) << " = "
//                    << int(
//                         repetitions.at(0) * repetitions.at(1) * repetitions.at(2)
//                         * pow(2, 3 * ref_level) * p * p * p) << endl
//                    << "y+ (wrt. cells): "  << yplus << "   dx+ = " << dxplus << ", "
//                    << "dz+ = " << dzplus << endl << "          -----         " << endl
//                    << "          -----         " << endl;
//                    << "==================================================="
//                        << endl << "               dt  = " << dt << endl << "               dt+ = "
//                        << dt/(viscosity/utau/utau) << endl << "    u_tau cross time = "
//                        << length / utau << " = "
//                        << int(length / utau / dt) << " steps" << endl
//                        << "===================================================" << endl
//                        << endl

    boost::shared_ptr<ProblemDescription<3> > mixingLayer =
            boost::make_shared<MixingLayer3D>(viscosity, refinement_level, meshname, randuscaling, randuname, U * uscaling);
    /////////////////////////////////////////////////
    // run solver
    //////////////////////////////////////////////////
    CompressibleCFDSolver<3> solver(configuration, mixingLayer);
    const size_t table_output_lines_per_10s = 300;
    configuration->setOutputTableInterval(1 + 10.0 / solver.getTimeStepSize() / table_output_lines_per_10s);
    solver.appendDataProcessor(boost::make_shared<ShearLayerStats>(solver, configuration->getOutputDirectory(), shearLayerThickness, Re));
    solver.run();
    pout << "step-mixingLayer terminated." << endl;
    return 0;
}
