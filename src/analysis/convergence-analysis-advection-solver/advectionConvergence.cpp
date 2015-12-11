/**
 * @file advectionConvergence.cpp
 * @short Checks the convergence of the SEDG advection solver, analyzes the local discretization error in dependence of dt, dx and fe_order
 * @date 13.02.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "natrium/benchmarks/AdvectionBenchmark.h"

using namespace natrium;
using namespace natrium::AdvectionBenchmark;

const bool OUTPUT = false;
const bool SMOOTH = true;

// Main function
// Test the dependence between dt,dx,p and the global discretization error
int main() {

	MPIGuard::getInstance();

	pout << "Starting convergence test for the SEDG advection solver.." << endl;

	// Make results dir
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/convergence-advection-solver";
	if (mkdir(dirName.str().c_str(), 0777) == -1) {
		if (errno != EEXIST) {
			cerr << "Fehler in mkdir: " << strerror(errno) << endl;
		}
	}

	std::stringstream resFileName;
	resFileName << dirName.str() << "/results.txt";
	std::ofstream out(resFileName.str().c_str());

	size_t numberOfTimeSteps = 0;
	out << "#Number of time steps: " << numberOfTimeSteps << endl;
	out << "#| level | " << "p  | " << "dx    | " << "dt     | "
			<< "error (2-norm)  | " << "error (sup-norm) | "
			<< " time/step (sec) |" << "non-0 in SparsityPattern" << endl;

	for (size_t refinementLevel = 1; refinementLevel < 5; refinementLevel++) {
		for (size_t feOrder = 1; feOrder < 12; feOrder++) {
			pout << "N = " << refinementLevel << ",     " << "p = " << feOrder
					<< ", " << endl;
			double deltaX = 1. / (pow(2, refinementLevel));
			pout << "dx = " << deltaX << ",    " << endl;
			//for (int i = -1; i < 5; i++) {
			//	double deltaT = deltaX * pow(0.5, i);
			double deltaT = 0.4 * pow(0.5, refinementLevel)
					/ ((feOrder + 1) * (feOrder + 1));
			pout << "dt = " << deltaT << "; " << endl;
			numberOfTimeSteps = 1.0 / deltaT;
			if (numberOfTimeSteps <= 5) {
				continue;
			}
			AdvectionResult result = oneTest(refinementLevel, feOrder, deltaT,
					numberOfTimeSteps, SMOOTH, true, false);
			out << result.refinementLevel << " " << result.fe_order << " "
					<< result.deltaX << " " << result.deltaT << " "
					<< result.norm2 << " " << result.normSup << " "
					<< result.timesec;
			//}
			//out << endl;
			pout << "----------" << endl;
		}
		pout << endl;
		out << endl;
	}

	pout << "done." << endl;
} /* Main */
