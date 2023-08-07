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
    parser.setArgument<int>("nout", "output vtk every nout steps", 1000);
    parser.setArgument<int>("nstats", "output stats every nstats steps", 20);
    parser.setArgument<int>("squash", "squash grid towards centre", 0);
    parser.setArgument<int>("print", "print calculations of initial velocity", 0);
    parser.setArgument<int>("recalculate", "recalculate initial velocity", 1);
    parser.setPositionalArgument<int>("ref-level",
                                      "Refinement level of the computation grid.");
    parser.setArgument<int>("grid-repetitions",
                            "Number of grid cells along each axis before global refinement; "
                            "to produce grids with refinements that are not powers of two.", 3); // TODO 1,2,3
    try { parser.importOptions();
    } catch (HelpMessageStop&) { return 0;
    }
    double Re = parser.getArgument<int>("Re");
    double refinement_level = parser.getArgument<int>("ref-level");
    long nout = parser.getArgument<int>("nout");
    auto time = parser.getArgument<double>("time");
    bool squash = parser.getArgument<int>("squash");
    bool print = parser.getArgument<int>("print");
    bool recalculate = parser.getArgument<int>("recalculate");

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
    const double cs = U / Ma;

    // chose scaling so that the right Ma-number is achieved
    const double scaling = sqrt(3) * cs;

    // setup configuration
    boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<SolverConfiguration>();
    configuration->setUserInteraction(false);
    configuration->setOutputCheckpointInterval(nout*100);
    configuration->setOutputSolutionInterval(nout);
    configuration->setSimulationEndTime(time);
    configuration->setOutputGlobalTurbulenceStatistics(false);
    configuration->setOutputCompressibleTurbulenceStatistics(false);
    configuration->setOutputShearLayerStatistics(false);
    configuration->setOutputShearLayerInterval(parser.getArgument<int>("nstats"));
    configuration->setStencilScaling(scaling);
    configuration->setStencil(Stencil_D3Q45);
    configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
    configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
    configuration->setHeatCapacityRatioGamma(1.4);
    configuration->setPrandtlNumber(0.71);
    configuration->setSedgOrderOfFiniteElement(3); // TODO: default, 2, 3
//    configuration->setInitializationScheme(COMPRESSIBLE_ITERATIVE);

    parser.applyToSolverConfiguration(*configuration);

    // standard output dir
    string m_dirname;
    if (not parser.hasArgument("output-dir")){
        std::stringstream dirName;
        dirName << getenv("NATRIUM_HOME") << "/step-mixingLayer/Re" << Re
                << "-Ma" << Ma
                << "-ref" << refinement_level
                << "-p" << configuration->getSedgOrderOfFiniteElement();
//        dirName << "-coll" << static_cast<int>(configuration->getCollisionScheme())
//                << "-sl" << static_cast<int>(configuration->getAdvectionScheme())
//                << "-";
        if (configuration->getAdvectionScheme() != SEMI_LAGRANGIAN)
            dirName << "-int" << static_cast<int>(configuration->getTimeIntegrator()) << "_" << static_cast<int>(configuration->getDealIntegrator());
        dirName << "-CFL" << configuration->getCFL();
//        dirName << "-sten" << static_cast<int>(configuration->getStencil());
//        if (configuration->isFiltering())
//            dirName << "-filt" << static_cast<int>(configuration->getFilteringScheme()) << "by_max_degree";
        if (configuration->getRegularizationScheme() != NO_REGULARIZATION)
            dirName << "-reg" << static_cast<int>(configuration->getRegularizationScheme());
//        if (configuration->getEquilibriumScheme()!= BGK_EQUILIBRIUM)
//            dirName << "-equili" << static_cast<int>(configuration->getEquilibriumScheme());
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

    boost::shared_ptr<ProblemDescription<3> > mixingLayer =
            boost::make_shared<MixingLayer3D>(viscosity, refinement_level, squash, print, recalculate, m_dirname, U);
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
