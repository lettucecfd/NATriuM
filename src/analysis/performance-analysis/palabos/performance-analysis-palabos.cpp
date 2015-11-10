#include <time.h>
//#include <ctime>

#include "palabos2D.h"
#include "palabos2D.hh"
#include <iostream>
#include <iomanip>
#include <math.h>

#include "../palabos/pramsmod.h"	// includes modified class for parameters

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor
#define _USE_MATH_DEFINES

//Global parameters. They will be redefined in relation to config.xml in advance of each computation
ModIncomprFlowParam<T> prams((T) 0,  	// uMax
		(T) 0,  		// Re
		0,       	// N
		0,        // lx
		0         // ly
		);

T rho0 = 1.; // All cells have initially density rho = 1

//Describes the analytical solution of the Taylor-Green-Vortex at a given time and location
void analyticalSolutionU(plint iX, plint iY, plint iteratedSteps,
		Array<T, 2>& u) {
	T factor = (-2) * (iteratedSteps) * prams.getDeltaT()
			* prams.getPhysicalNu();

	u[0] = sin(2 * M_PI * iX / prams.getResolution())
			* cos(2 * M_PI * iY / prams.getResolution()) * exp(factor);
	u[1] = -cos(2 * M_PI * iX / prams.getResolution())
			* sin(2 * M_PI * iY / prams.getResolution()) * exp(factor);
}

//Initializes the start velocities of the Taylor-Green-Vortex by using the analytical solution
void velocityOfVortex(plint iX, plint iY, T& rho, Array<T, 2>& u) {

	analyticalSolutionU(iX, iY, 0, u); // Using the analytical solution as the starting parameters for iT=0

	u[0] = u[0] * prams.getLatticeU();
	u[1] = u[1] * prams.getLatticeU();
	rho = rho0;

}

void initializeLattice(MultiBlockLattice2D<T, DESCRIPTOR>& lattice)
{
	initializeAtEquilibrium (
	lattice, lattice.getBoundingBox(), velocityOfVortex );

	lattice.initialize();
}

int main(int argc, char* argv[]) {

	string outputFolder, resultFile;
	plint n_min, n_max, k_min, k_max;
	T max_time, u_max, reynolds, physicalSize;

	try {
		// Open the XMl configuration file.
		XMLreader xmlFile("config.xml");
		xmlFile["Data"]["outputFolder"].read(outputFolder);
		xmlFile["Data"]["resultFile"].read(resultFile);

		xmlFile["Geometry"]["size"]["n_min"].read(n_min);
		xmlFile["Geometry"]["size"]["n_max"].read(n_max);
		xmlFile["Geometry"]["size"]["k_min"].read(k_min);
		xmlFile["Geometry"]["size"]["k_max"].read(k_max);
		xmlFile["Geometry"]["size"]["physicalSize"].read(physicalSize);

		//   xmlFile["Geometry"]["viscosity"].read(nu_lat);
		xmlFile["Geometry"]["reynolds"].read(reynolds);
		xmlFile["Geometry"]["u_max"].read(u_max);
		xmlFile["Geometry"]["time"].read(max_time);
	}

	catch (PlbIOException& exception) {
		pcout << exception.what() << endl;
		return -1;

	}

	plbInit(&argc, &argv);

	global::directories().setOutputDir(outputFolder);

	std::stringstream filename_plb;
	time_t now = time(0);
	tm *ltm = localtime(&now);
	filename_plb << getenv("NATRIUM_HOME")
			<< "/performance-analysis/performance-Plb-" << getenv("HOSTNAME")
			<< "-" << 1900 + ltm->tm_year << "-" << ltm->tm_mon + 1 << "-"
			<< ltm->tm_mday << "_" << ltm->tm_hour << "-" << ltm->tm_min
			<< ".txt";
	plb_ofstream outfile_plb(filename_plb.str().c_str());
	outfile_plb
			<< "# refinement   p    #dofs   #steps    init time (sec)   loop time (sec)   time for one iteration (sec)"
			<< endl;

	prams.Re = reynolds;
	prams.lx = physicalSize;
	prams.ly = physicalSize;
	//prams.physicalSize = physicalSize;

	//Compute in relation to k*2^n
	for (int k = k_min; k <= k_max; k++) {
		for (int n = n_min; n <= n_max; n++) {
			cout << "with N = " << n << "; p = " << k << " ... " << endl;

			double time1, time2, timestart;

			//// INIT
			timestart = clock();

			plint n_lat = (k + 1) * pow(2, n); // Grid size is computed by k*2^n
			prams.resolution = n_lat;
			prams.latticeU = u_max;
			MultiBlockLattice2D<T, DESCRIPTOR> lattice (n_lat, n_lat, new BGKdynamics<T,DESCRIPTOR>(prams.getOmega()) );
			lattice.periodicity().toggleAll(true); // Set periodic boundaries.
			lattice.toggleInternalStatistics(false); // no output
			initializeLattice(lattice); // Initializes the lattice using the starting parameters
			time1 = clock() - timestart;

			// LOOP

			plint steps = 0;
			plint n_steps = 10000 / pow(2,n);

			for (steps = 0; steps < n_steps; steps++) // Loop until max_time is reached
					{

				// Execute lattice Boltzmann iteration.
				lattice.collideAndStream();

			}
			time2 = clock() - time1 - timestart;
			time1 /= CLOCKS_PER_SEC;
			time2 /= CLOCKS_PER_SEC;

			outfile_plb << setprecision(10) << n << " " << k << " "
					<< pow((k + 1) * pow(2, n), 2) << " " << steps << " "
					<< time1 << " " << time2 << " "
					<< time2 / n_steps << endl;

		}
		outfile_plb << endl;

	}

}

