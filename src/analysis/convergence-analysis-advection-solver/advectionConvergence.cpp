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

#include "deal.II/numerics/data_out.h"

#include "advection/SEDGMinLee.h"

#include "problemdescription/ProblemDescription.h"

#include "timeintegration/RungeKutta5LowStorage.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "utilities/BasicNames.h"
#include "utilities/Math.h"

#include "PeriodicTestDomain2D.h"

using namespace natrium;

//#define OUTPUT

// analytic solution (where should the bump be after a certain time)
void getAnalyticSolution(double time, distributed_vector& analyticSolution,
		const vector<dealii::Point<2> >& supportPoints) {
	assert(analyticSolution.size() == supportPoints.size());
	assert(supportPoints.size() > 0);

	//dealii::Point<2> midPoint(0.25, 0.25);
	//double R = 0.125;
	double lambda = 1;

	dealii::Point<2> originalPoint;
	for (size_t i = 0; i < supportPoints.size(); i++) {
		originalPoint = supportPoints.at(i);
		// move back to original point
		originalPoint(0) -= time;
		// transform back to [0,1]
		while (originalPoint(0) < 0.0) {
			originalPoint(0) += 1.0;
		}
		while (originalPoint(0) > 1.0) {
			originalPoint(0) -= 1.0;
		}
		/*double distance = originalPoint.distance(midPoint);
		 if (distance <= 2*R) {
		 // smooth function!
		 double smooth = 0.5 - 1.0/(4*atan(1)) * atan(1.0/pow(5*(R-distance),2) - 1.0/pow(5*distance,2));
		 //double smooth = 0.1*cos(Math::PI / 0.5 * distance);
		 analyticSolution(i) = 1 + smooth;//; + 0.1;/smooth;//cos(Math::PI / 0.5 * distance);  1 + 0.1;
		 } else {
		 analyticSolution(i) = 1;
		 }*/
		double h = sin(8 * atan(1) * originalPoint(0) * lambda);
		analyticSolution(i) = ((h > 0) - (h < 0)) * pow(fabs(h), 2.0);
		//analyticSolution(i) = exp(-(originalPoint(0)-0.25)*(originalPoint(0)-0.25));
	}
}

std::string oneTest(size_t refinementLevel, size_t fe_order, double deltaT,
		size_t numberOfTimeSteps, bool useCentralFlux = false) {

	double deltaX = 1. / (pow(2, refinementLevel));

	// create problem and solver
	PeriodicTestDomain2D periodic(refinementLevel);
	SEDGMinLee<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), "", useCentralFlux);
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();

	// create smooth initial conditions
	distributed_vector f(streaming.getNumberOfDoFs());
	// zero-vector for time integrator
	distributed_vector g(streaming.getNumberOfDoFs());
	vector<dealii::Point<2> > supportPoints(
			streaming.getDoFHandler()->n_dofs());
	dealii::DoFTools::map_dofs_to_support_points(streaming.getMapping(),
			*streaming.getDoFHandler(), supportPoints);
	getAnalyticSolution(0.0, f, supportPoints);
	assert(g.all_zero());

	// RK5 for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern(0));
	advectionMatrix.copy_from(matrices.block(0, 0));
	RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector> RK5(
			deltaT, f.size());

	distributed_vector fAnalytic(f.size());

	double timesec, timestart;
	timestart = clock();

#ifdef OUTPUT
	// prepare output directory
	std::stringstream dirName;
	dirName << getenv("NATRIUM_HOME") << "/convergence-advection-solver/Level_"
	<< refinementLevel << "_p_" << fe_order;
	if (mkdir(dirName.str().c_str(), 0777) == -1) {
		if (errno != EEXIST) {
			cerr << "Fehler in mkdir: " << strerror(errno) << endl;
		}
	}

	for (size_t i = 0; i < numberOfTimeSteps; i++) {
		//stream
		RK5.step(f, advectionMatrix, g);

		// output
		if (i % 10 == 0) {
			std::stringstream str;
			str << dirName.str().c_str() << "/t_"
			<< i << ".vtu";
			std::string filename = str.str();
			std::ofstream vtu_output(filename.c_str());
			dealii::DataOut<2> data_out;
			data_out.attach_dof_handler(*streaming.getDoFHandler());
			data_out.add_data_vector(f, "f");
			getAnalyticSolution(deltaT * i, fAnalytic, supportPoints);
			data_out.add_data_vector(fAnalytic, "f_ref");
			data_out.build_patches(20);
			data_out.write_vtu(vtu_output);
		}

	}
#else

	for (size_t i = 0; i < numberOfTimeSteps; i++) {
		//stream
		RK5.step(f, advectionMatrix, g);
	}

#endif
	timesec = clock() - timestart;
	timesec /= numberOfTimeSteps;
	timesec /= CLOCKS_PER_SEC;

	// compare with analytic solution by sup-norm and 2-norm
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
	result << "      " << refinementLevel << "      " << fe_order << "     "
			<< deltaX << "     " << deltaT << "      " << norm2 << "       "
			<< normSup << "    " << timesec << "    "
			<< advectionMatrix.get_sparsity_pattern().n_nonzero_elements();

	std::stringstream fileName;
	fileName << getenv("NATRIUM_HOME") << "/convergence-advection-solver/Level_"
			<< refinementLevel << "_p_" << fe_order << ".sp";
	std::ofstream sp_file(fileName.str().c_str());
	streaming.getBlockSparsityPattern().print_gnuplot(sp_file);
	return result.str();

}

// Main function
// Test the dependence between dt,dx,p and the global discretization error
int main() {
	cout << "Starting convergence test for the SEDG advection solver.." << endl;

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
			cout << "N = " << refinementLevel << ",     " << "p = " << feOrder
					<< ", " << endl;
			double deltaX = 1. / (pow(2, refinementLevel));
			cout << "dx = " << deltaX << ",    " << endl;
			//for (int i = -1; i < 5; i++) {
			//	double deltaT = deltaX * pow(0.5, i);
			double deltaT = 0.4 * pow(0.5, refinementLevel)
					/ (feOrder * feOrder);
			cout << "dt = " << deltaT << "; " << endl;
			numberOfTimeSteps = 1.0 / deltaT;
			if (numberOfTimeSteps <= 5) {
				continue;
			}
			std::string result = oneTest(refinementLevel, feOrder, deltaT,
					numberOfTimeSteps);
			out << result.c_str() << endl;
			//}
			//out << endl;
			cout << "----------" << endl;
		}
		cout << endl;
		out << endl;
	}

	cout << "done." << endl;
} /* Main */
