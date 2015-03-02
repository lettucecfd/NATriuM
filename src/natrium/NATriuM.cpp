//============================================================================
// Name        : NATriuM.cpp
// Author      : Andreas Kraemer, andreas.kraemer@h-brs.de
// Version     : 0
// Copyright   : University of Applied Sciences Sankt Augustin, Germany
// Description : Numerics and Algorithms for Tribology
//============================================================================


/*! \mainpage NATriuM - A discontinuous Galerkin Lattice Boltzmann simulation tool based on deal.II
 * \section lbm_sec The Lattice Boltzmann Method
 *
 * Lattice Boltzmann (LB) methods are an alternative approach for the simulation of fluid
 * flows. Modeling particulate processes (advection and collision) they indirectly solve the
 * fluid mechanical Navier-Stokes equations. Since 1988 LB has established besides classical
 * finite element, finite volume and finite difference schemes, particularly for the simulation
 * of complex flows. Concretely, they provide a straightforwardmodeling of multiphase flows
 * and flows in complex geometries. Both kinds of simulations are nowadays indispensable,
 * especially in engineering applications. The LB method has proven suitable for simulating
 * complex flows on parallel computers.
 * However, the original LB algorithm is limited to regular computation grids (the lattices),
 * which is a severe shortcoming, specifically in resolving different scales and flows at high
 * Reynolds numbers (Re). Furthermore, it is restricted to at most second order accuracy in time and space
 * due to the inherent Euler time integration. Consequently, researchers attempted to transfer the LB concept
 * to irregular grids and higher order accuracy. Various approaches have been published.
 *
 * \section sedg_sec LBM on unstructured grids
 * A very promising approach for LBM on irregular grids was published by Lee and Lin in 2003. They managed
 * to split the Boltzmann equation into collision and advection, without restricting it to regular grids.
 * They could solve the advection equation (streaming step) using a Finite Element scheme on a triangular
 * grid with high order and even implicit time integrators. Their method was picked up by Min and Lee in
 * 2011, who used a spectral element discontinuous Galerkin (SEDG) solver for the advection and obtained
 * a high-order convergent SEDG-LBM scheme. NATriuM is based on this approach. Obviously, the sophisticated
 * streaming procedure comes at the price of higher computation time. We have measured factors of ~60 in
 * a comparison with the Open Source LBM software Palabos (for BGK-collision and equal CFL numbers).
 *
 * \section tpl_sec Third party libraries (TPLs)
 * NATriuM inherits a huge part of its functionality from third party libraries, most importantly deal.II.
 * Deal.II is a finite element library which was awarded the Wilkinson price for numerical software in 2007.
 * Together with a huge FEM functionality, it offers wrappers to other TPLs like Trilinos to speed
 * up simulations on multicore architectures. We use the Trilinos data structures (sparse matrices and vectors).
 * The boost unit test module is used to structure the testing code. CMake is used for user-friendly compilation.
 *
 * \section struct_sec Code structure
 * At the moment, NATriuM can only be used as a library, but it is going to be usable from the command line
 * in a later stage of the project. A class diagram is in the documentation folder, but it could be slightly
 * outdated. The source directory contains the following folders:
 *   - natrium: NATriuM library. Contains the code for the SEDG-Solver.
 *   - test: Unit and integration tests.
 *   - doc: Code documentation. Most importantly, this Doxygen documentation.
 *   - examples: Example simulations that use the NATriuM code.
 *   - analysis: Scripts for numerical analysis of the code (convergence tests, etc.)
 *   - preprocessing: Parameter GUI.
 *   - postprocessing: Nothing special (at the moment).
 *
 * \section workflow_sec Workflow
 * It is recommended to use the code in connection with sophisticated pre- and postprocessing tools.
 *   - Preprocessing: Deal.II supports various mesh grid formats. E.g. Salome can be used for grid generation.
 *   - Postprocessing: Paraview is highly recommended.
 *
 * The Open Fuel Cell simulator (OpenFCST) also uses the workflow Salome->deal.II->Paraview. It might be interesting
 * for you to take a look at the OpenFCST user manual (online available) which includes for example a short Tutorial
 * on interfacing Salome and deal.II.
 *
 * \section start_sec Getting started
 * For simple applications, the triangulation can be generated within deal.II. The examples in the Examples section will help
 * you to get grip on the code.
 *
 * \section install_sec Installation
 *
 * The installation of NATRiuM requires several third party libraries. For a detailed documentation see the file INSTALL.txt.
 *
 */


#include "deal.II/grid/tria.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "solver/CFDSolver.h"

#include "utilities/BasicNames.h"

using namespace natrium;
//using dealii::Triangulation;


int main() {

	std::ifstream iFile("/home/kraemer/test.txt");

	cout << "Starting NATriuM." << endl;

	std::string a;
	cout << "gimme: ";
	do {
	iFile >> a;
	//Triangulation<2> tri;

	std::cout << a << endl;
	} while (not iFile.eof());
	vector<double> v(2,0.0);
	D2Q9IncompressibleModel boltzmannModel;

	//CFDSolver<2>* solver = new CFDSolver<2>();



	//double rho = 1;
	//cout << boltzmannModel->getEquilibriumDistribution(0,v,rho) << endl;

	cout << "NATriuM terminated." << endl;

	return 0;

}
