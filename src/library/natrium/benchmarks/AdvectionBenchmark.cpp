/*
 * AdvectionBenchmark.cpp
 *
 *  Created on: May 11, 2015
 *      Author: kraemer
 */

#include "AdvectionBenchmark.h"

#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "deal.II/numerics/data_out.h"

#include "../advection/SEDGMinLee.h"

#include "../problemdescription/ProblemDescription.h"

#include "../timeintegration/RungeKutta5LowStorage.h"

#include "../stencils/D2Q9.h"

#include "../utilities/Math.h"

#include "../benchmarks/PeriodicTestDomain2D.h"


namespace natrium {

template class SEDGMinLee<2>;
template class RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector>;

namespace AdvectionBenchmark {
// analytic solution (where should the bump be after a certain time)
void getAnalyticSolution(double time, distributed_vector& analyticSolution,
		const vector<dealii::Point<2> >& supportPoints, bool is_smooth) {
	assert(analyticSolution.size() == supportPoints.size());
	assert(supportPoints.size() > 0);

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
		if (is_smooth) {
			// see Hesthaven, p. 27
			double h = sin(8 * atan(1) * originalPoint(0) * lambda);
			analyticSolution(i) = h;
		} else {
			// see Hesthaven, p. 83
			double h = sin(8 * atan(1) * originalPoint(0) * lambda);
			analyticSolution(i) = ((h > 0) - (h < 0)) * pow(fabs(h), 2.0);
		}
	} /* for i in points */
} /* getAnalyticSolution */

AdvectionResult oneTest(size_t refinementLevel, size_t fe_order, double deltaT,
		size_t numberOfTimeSteps, bool is_smooth, bool output_to_std_dir,
		bool useCentralFlux) {

	AdvectionResult result;
	result.refinementLevel = refinementLevel;
	result.fe_order = fe_order;
	result.deltaT = deltaT;

	result.deltaX = 1. / (pow(2, refinementLevel));

	// create problem and solver
	PeriodicTestDomain2D periodic(refinementLevel);
	SEDGMinLee<2> streaming(periodic.getMesh(),
			periodic.getBoundaries(), fe_order, make_shared<D2Q9>(), "",
			useCentralFlux);
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();

#ifdef WITH_TRILINOS_MPI
	// create smooth initial conditions
	distributed_vector f(streaming.getLocallyOwnedDofs(), streaming.getLocallyRelevantDofs(), MPI_COMM_WORLD);
	// zero-vector for time integrator
	distributed_vector g(streaming.getLocallyOwnedDofs(), streaming.getLocallyRelevantDofs(), MPI_COMM_WORLD);
#else
	// create smooth initial conditions
	distributed_vector f(streaming.getNumberOfDoFs());
	// zero-vector for time integrator
	distributed_vector g(streaming.getNumberOfDoFs());
#endif
	vector<dealii::Point<2> > supportPoints(
			streaming.getDoFHandler()->n_dofs());
	dealii::DoFTools::map_dofs_to_support_points(streaming.getMapping(),
			*streaming.getDoFHandler(), supportPoints);
	getAnalyticSolution(0.0, f, supportPoints, is_smooth);
	assert(g.all_zero());

	// RK5 for streaming in direction (1,0)
#ifdef WITH_TRILINOS
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(matrices.block(0, 0));
	advectionMatrix.copy_from(matrices.block(0, 0));
#else
	// TODO (AK) Can we remove the ifdef?
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern(0));
	advectionMatrix.copy_from(matrices.block(0, 0));
#endif
	shared_ptr<TimeIntegrator<distributed_sparse_matrix,
	distributed_vector> > RK5 = make_shared<RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector> >(
			deltaT, f);

#ifdef WITH_TRILINOS_MPI
	distributed_vector fAnalytic(streaming.getLocallyOwnedDofs(), streaming.getLocallyRelevantDofs(), MPI_COMM_WORLD);
#else
	distributed_vector fAnalytic(f.size());
#endif

	double timestart;
	timestart = clock();

	if (output_to_std_dir) {
		// prepare output directory
		std::stringstream dirName;
		dirName << getenv("NATRIUM_HOME")
				<< "/convergence-advection-solver/Level_" << refinementLevel
				<< "_p_" << fe_order;
		if (mkdir(dirName.str().c_str(), 0777) == -1) {
			if (errno != EEXIST) {
				cerr << "Fehler in mkdir: " << strerror(errno) << endl;
			}
		}

		for (size_t i = 0; i < numberOfTimeSteps; i++) {
			//stream
			RK5->step(f, advectionMatrix, g);

			// output
			if (i % 10 == 0) {
				std::stringstream str;
				str << dirName.str().c_str() << "/t_" << i << ".vtu";
				std::string filename = str.str();
				std::ofstream vtu_output(filename.c_str());
				dealii::DataOut<2> data_out;
				data_out.attach_dof_handler(*streaming.getDoFHandler());
				data_out.add_data_vector(f, "f");
				getAnalyticSolution(deltaT * i, fAnalytic, supportPoints, is_smooth);
				data_out.add_data_vector(fAnalytic, "f_ref");
				data_out.build_patches(20);
				data_out.write_vtu(vtu_output);
			}

		}
	} else {

		for (size_t i = 0; i < numberOfTimeSteps; i++) {
			//stream
			RK5->step(f, advectionMatrix, g);
		}

	}
	result.timesec = clock() - timestart;
	result.timesec /= numberOfTimeSteps;
	result.timesec /= CLOCKS_PER_SEC;

	// compare with analytic solution by sup-norm and 2-norm
	getAnalyticSolution(deltaT * numberOfTimeSteps, fAnalytic, supportPoints, is_smooth);
	result.normSup = 0.0;
	result.norm2 = 0.0;
	for (size_t i = 0; i < f.size(); i++) {
		double error = fabs(f(i) - fAnalytic(i));
		if (error > result.normSup) {
			result.normSup = error;
		}
		result.norm2 += (error * error);
	}
	result.norm2 = sqrt(result.norm2);

	std::stringstream fileName;
	fileName << getenv("NATRIUM_HOME") << "/convergence-advection-solver/Level_"
			<< refinementLevel << "_p_" << fe_order << ".sp";
	std::ofstream sp_file(fileName.str().c_str());
	streaming.getBlockSparsityPattern().print_gnuplot(sp_file);
	return result;

} /* oneTest */

} /*namespace  AdvectionBenchmark */
} /* namespace natrium */
