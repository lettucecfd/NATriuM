/*
 * step-6.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: bajat
 */

#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/BenchmarkCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"

#include "natrium/stencils/D3Q19.h"

#include "natrium/problemdescription/Benchmark.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/CFDSolverUtilities.h"

#include "natrium/benchmarks/TaylorGreenVortex3D.h"

using namespace natrium;

// Main function
int main() {

  MPIGuard::getInstance();

  pout << "Starting NATriuM step-16 ..." << endl;

  /////////////////////////////////////////////////
  // set parameters, set up configuration object
  //////////////////////////////////////////////////

  // Re = viscosity/(2*pi)
  const double viscosity = 1;
  // C-E-approach: constant stencil scaling
  // specify Mach number
  const double Ma = 0.1;
  // zunaechst: fixed order of FE
  const double orderOfFiniteElement = 6;

  // chose scaling so that the right Ma-number is achieved
  double scaling = sqrt(3) * 1 / Ma;

  const double refinementLevel = 1;


  boost::shared_ptr<Benchmark<3> >  taylorGreen = boost::make_shared<
     TaylorGreenVortex3D>(viscosity, refinementLevel);


    // the scaling has to be orders of magnitude greater than the boundary velocity
    double dt = CFDSolverUtilities::calculateTimestep<3>(
					*taylorGreen->getMesh(), orderOfFiniteElement,
					D3Q19(scaling), 0.4);

    pout << "dt = " << dt << " ...";

    // time measurement variables
    double time1, time2, timestart;

    // setup configuration
    std::stringstream dirName;
    dirName << getenv("NATRIUM_HOME") << "/step-16";
    boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
        SolverConfiguration>();
    //configuration->setSwitchOutputOff(true);
    configuration->setOutputDirectory(dirName.str());
    configuration->setRestartAtLastCheckpoint(false);
    configuration->setUserInteraction(false);
    configuration->setOutputTableInterval(10);
    configuration->setOutputCheckpointInterval(1000);
    configuration->setOutputSolutionInterval(10);
    configuration->setNumberOfTimeSteps(10000000);
    configuration->setInitializationScheme(EQUILIBRIUM);
    configuration->setSedgOrderOfFiniteElement(orderOfFiniteElement);
    configuration->setStencilScaling(scaling);
    configuration->setStencil(Stencil_D3Q19);
    //configuration->setCommandLineVerbosity(BASIC);
    configuration->setTimeStepSize(dt);
    if (dt > 0.1) {
      pout << "Timestep too big." << endl;
    }

    //configuration->setNumberOfTimeSteps(1.0 / dt);

    timestart = clock();
    BenchmarkCFDSolver<3> solver(configuration, taylorGreen);
    time1 = clock() - timestart;


      solver.run();
      time2 = clock() - time1 - timestart;
      time1 /= CLOCKS_PER_SEC;
      time2 /= CLOCKS_PER_SEC;
      pout << " OK ... Init: " << time1 << " sec; Run: " << time2
          << " sec." << endl;



  pout << "step-1 terminated." << endl;

  return 0;

}
