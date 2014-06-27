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

#include "advection/SEDGMinLee.h"

#include "problemdescription/ProblemDescription.h"

#include "timeintegration/RungeKutta5LowStorage.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "utilities/BasicNames.h"
#include "utilities/Math.h"

#include "PeriodicTestDomain2D.h"

using namespace natrium;

// analytic solution (where should the bump be after a certain time)
void getAnalyticSolution(double time, distributed_vector& analyticSolution,
		const vector<dealii::Point<2> >& supportPoints) {
	assert(analyticSolution.size() == supportPoints.size());
	assert(supportPoints.size() > 0);

	dealii::Point<2> midPoint(0.25, 0.25);
	dealii::Point<2> originalPoint;
	for (size_t i = 0; i < supportPoints.size(); i++) {
		originalPoint = supportPoints.at(i);
		// move back to original point
		originalPoint(0) -= time;
		// transform back to [0,1]
		while (originalPoint(0) < 0.0) {
			originalPoint(0) += 1.0;
		}
		while (originalPoint(1) > 1.0) {
			originalPoint(0) -= 1.0;
		}
		double distance = originalPoint.distance(midPoint);
		if (distance <= 0.25) {
			analyticSolution(i) = 1 + 0.1; //* cos(Math::PI / 0.5 * distance);
		} else {
			analyticSolution(i) = 1;
		}
	}
}

std::string oneTest(size_t refinementLevel, size_t fe_order, double deltaT,
		size_t numberOfTimeSteps, bool useCentralFlux = false) {

	double deltaX = 1. / (pow(2, refinementLevel) * (fe_order - 1));

	// create problem and solver
	PeriodicTestDomain2D periodic(refinementLevel);
	SEDGMinLee<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), "", useCentralFlux);
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();

	// create smooth initial conditions
	distributed_vector f(streaming.getNumberOfDoFs());
	vector<dealii::Point<2> > supportPoints(
			streaming.getDoFHandler()->n_dofs());
	dealii::DoFTools::map_dofs_to_support_points(streaming.getMapping(),
			*streaming.getDoFHandler(), supportPoints);
	getAnalyticSolution(0.0, f, supportPoints);

	// RK5 for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern(0));
	advectionMatrix.copy_from(matrices.block(0,0));
	RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector> RK5(deltaT, f.size());

	for (size_t i = 0; i < numberOfTimeSteps; i++) {
		//stream
		RK5.step(f, advectionMatrix);
	}

	// compare with analytic solution by sup-norm and 2-norm
	distributed_vector fAnalytic(f.size());
	getAnalyticSolution(deltaT * numberOfTimeSteps, fAnalytic, supportPoints);
	double normSup = 0.0;
	double norm2 = 0.0;
	for (size_t i = 0; i < f.size(); i++) {
		double error = fabs(f(i) - fAnalytic(i));
		if (error > normSup) {
			normSup = error;
		}
		norm2 += (error * error);
	}
	norm2 = sqrt(norm2);

	std::stringstream result;
	result << "      " << refinementLevel << "      " << fe_order << "  "
			<< deltaX << " " << deltaT << "      " << norm2 << "       "
			<< normSup;
	return result.str();

}

// Main function
// Test the dependence between dt,dx,p and the global discretization error
int main() {
	cout << "Starting convergence test for the SEDG advection solver.."
			<< endl;

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
			<< "error (2-norm)  | " << "error (sup-norm) | " << endl;


	for (size_t feOrder = 2; feOrder < 6; feOrder++) {
		cout << "p = " << feOrder << ", " << endl;
		for (size_t refinementLevel = 3; refinementLevel < 7;
				refinementLevel++) {
			double deltaX = 1. / (pow(2, refinementLevel) * (feOrder - 1));
			cout << endl << "dx = " << deltaX << ",    " << endl;
			for (int i = -1; i < 5; i++) {
				double deltaT = deltaX * pow(0.5, i);
				cout << "dt = " << deltaT << "; " << endl;
				numberOfTimeSteps = 0.5/deltaT;
				std::string result = oneTest(refinementLevel, feOrder, deltaT,
						numberOfTimeSteps);
				out << result.c_str() << endl;
			}
			out << endl;
		}
		cout << endl;
	}

	cout << "done." << endl;
} /* Main */