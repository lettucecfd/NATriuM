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
#include "../advection/SemiLagrangian.h"

#include "../problemdescription/ProblemDescription.h"

#include "../timeintegration/ThetaMethod.h"
#include "../timeintegration/RungeKutta5LowStorage.h"
#include "../timeintegration/ExponentialTimeIntegrator.h"
#include "../timeintegration/DealIIWrapper.h"

#include "../stencils/D2Q9.h"

#include "../utilities/Math.h"

#include "../benchmarks/PeriodicTestDomain2D.h"

namespace natrium {

template<size_t dim> class SEDGMinLee;
template<size_t dim> class SemiLagrangian;

boost::shared_ptr<TimeIntegrator<distributed_sparse_matrix, distributed_vector> > make_integrator(
		TimeIntegratorName integrator_name,
		DealIntegratorName deal_integrator_name, double delta_t,
		distributed_vector& f) {

	/// Build time integrator
	switch (integrator_name) {
	case RUNGE_KUTTA_5STAGE: {
		return boost::make_shared<
				RungeKutta5LowStorage<distributed_sparse_matrix,
						distributed_vector> >(delta_t, f);
	}
	case THETA_METHOD: {
		const double theta = 0.5;
		return boost::make_shared<
				ThetaMethod<distributed_sparse_matrix, distributed_vector> >(
				delta_t, f, theta);
	}
	case EXPONENTIAL: {
		return boost::make_shared<
				ExponentialTimeIntegrator<distributed_sparse_matrix,
						distributed_vector> >(delta_t);
	}
	case OTHER: {
		// create configuration object to get standard preferences
		SolverConfiguration config;
		if (deal_integrator_name < 7) {
			return boost::make_shared<
					DealIIWrapper<distributed_sparse_matrix, distributed_vector> >(
					delta_t, deal_integrator_name, config.getDealLinearSolver());
		} else if (deal_integrator_name < 12) {
			return boost::make_shared<
					DealIIWrapper<distributed_sparse_matrix, distributed_vector> >(
					delta_t, deal_integrator_name, config.getDealLinearSolver(),
					config.getEmbeddedDealIntegratorCoarsenParameter(),
					config.getEmbeddedDealIntegratorRefinementParameter(),
					0.1 * delta_t, // minimum time step
					10 * delta_t, // maximum time step
					config.getEmbeddedDealIntegratorRefinementTolerance(),
					config.getEmbeddedDealIntegratorCoarsenTolerance());
		}
		assert(false);
		return NULL;
	}/* case OTHER */
	default: {
		assert(false);
		return NULL;
	}
	} /* switch integrator_name */
}

namespace AdvectionBenchmark {

// analytic solution (where should the bump be after a certain time)
double analytic_solution(double time, const dealii::Point<2>& x,
		bool is_smooth) {

	double lambda = 1;
	dealii::Point<2> originalPoint;
	originalPoint = x;
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
		return sin(8 * atan(1) * originalPoint(0) * lambda);
	} else {
		// see Hesthaven, p. 83
		double h = sin(8 * atan(1) * originalPoint(0) * lambda);
		return ((h > 0) - (h < 0)) * pow(fabs(h), 2.0);
	}
} /* analytic_solution */

void getAnalyticSolution(double time, distributed_vector& analyticSolution,
		const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints,
		const AdvectionOperator<2>& streaming, bool is_smooth) {
	// get Function instance
	const unsigned int dofs_per_cell = streaming.getFe()->dofs_per_cell;
	vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			streaming.getDoFHandler()->begin_active(), endc =
			streaming.getDoFHandler()->end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			cell->get_dof_indices(local_dof_indices);
			for (size_t i = 0; i < dofs_per_cell; i++) {
				assert(
						analyticSolution.in_local_range(
								local_dof_indices.at(i)));
				assert(
						supportPoints.find(local_dof_indices.at(i))
								!= supportPoints.end());
				analyticSolution(local_dof_indices.at(i)) = analytic_solution(
						time, supportPoints.at(local_dof_indices.at(i)),
						is_smooth);
			}
		} /* if is locally owned */
	} /* for all cells */
}

AdvectionResult oneTest(size_t refinementLevel, size_t fe_order, double deltaT,
		double t_end, const TimeIntegratorName integrator,
		const DealIntegratorName deal_integrator, bool is_smooth,
		bool semi_lagrangian, bool output_to_std_dir, bool useCentralFlux) {

	AdvectionResult result;
	result.refinementLevel = refinementLevel;
	result.fe_order = fe_order;
	result.deltaT = deltaT;

	result.deltaX = 1. / (pow(2, refinementLevel));

	// create problem and solver
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();

	boost::shared_ptr<AdvectionOperator<2> > s;
	if (semi_lagrangian) {
		s = boost::make_shared<SemiLagrangian<2> >(periodic.getMesh(),
				periodic.getBoundaries(), fe_order, boost::make_shared<D2Q9>(),
				deltaT);
	} else {
		s = boost::make_shared<SEDGMinLee<2> >(periodic.getMesh(),
				periodic.getBoundaries(), fe_order, boost::make_shared<D2Q9>(),
				useCentralFlux);
	}
	AdvectionOperator<2> & streaming = *s;
	streaming.setupDoFs();
	streaming.reassemble();
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();

#ifdef WITH_TRILINOS_MPI
	// create smooth initial conditions
	distributed_vector f(streaming.getLocallyOwnedDofs(), MPI_COMM_WORLD);
	// zero-vector for time integrator / tmp-vector for semi-lagrangian
	distributed_vector g(streaming.getLocallyOwnedDofs(), MPI_COMM_WORLD);
#else
	// create smooth initial conditions
	distributed_vector f(streaming.getNumberOfDoFs());
	// zero-vector for time integrator
	distributed_vector g(streaming.getNumberOfDoFs());
#endif
	map<dealii::types::global_dof_index, dealii::Point<2> > supportPoints;
	streaming.mapDoFsToSupportPoints(supportPoints);
	getAnalyticSolution(0.0, f, supportPoints, streaming, is_smooth);
	assert(g.all_zero());

	// RK5 for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(matrices.block(0, 0));
	advectionMatrix.copy_from(matrices.block(0, 0));
	boost::shared_ptr<
			TimeIntegrator<distributed_sparse_matrix, distributed_vector> > time_stepper;
	time_stepper = make_integrator(integrator, deal_integrator, deltaT, f);

	distributed_vector fAnalytic(streaming.getLocallyOwnedDofs(),
	MPI_COMM_WORLD);

	double timestart;
	timestart = clock();
	double t = 0;
	size_t i = 0;

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

		while (t < t_end) {
			i++;
			//for (size_t i = 0; i < numberOfTimeSteps; i++) {
			//stream
			if (semi_lagrangian) {
				g = f;
				t += deltaT;
				advectionMatrix.vmult(f, g);
			} else {
				t += time_stepper->step(f, advectionMatrix, g, t,
						time_stepper->getTimeStepSize());
			}
			// output
			if (i % 10 == 0) {
				std::stringstream str;
				str << dirName.str().c_str() << "/t_" << i << ".vtu";
				std::string filename = str.str();
				std::ofstream vtu_output(filename.c_str());
				dealii::DataOut<2> data_out;
				data_out.attach_dof_handler(*streaming.getDoFHandler());
				data_out.add_data_vector(f, "f");
				getAnalyticSolution(deltaT * i, fAnalytic, supportPoints,
						streaming, is_smooth);
				data_out.add_data_vector(fAnalytic, "f_ref");
				data_out.build_patches(20);
				data_out.write_vtu(vtu_output);
			}

		}
	} else {

		while (t < t_end) {
			i++;

			//stream
			if (semi_lagrangian) {
				g = f;
				t += deltaT;
				advectionMatrix.vmult(f, g);
			} else {
				t = time_stepper->step(f, advectionMatrix, g, t,
						time_stepper->getTimeStepSize());
			}
			//pout << t << endl;
		}

	}
	result.timesec = clock() - timestart;
	result.timesec /= i;
	result.timesec /= CLOCKS_PER_SEC;

	// compare with analytic solution by sup-norm and 2-norm
	getAnalyticSolution(t, fAnalytic, supportPoints, streaming, is_smooth);
	fAnalytic.add(-1, f);
	result.normSup = fAnalytic.linfty_norm();
	result.norm2 = fAnalytic.l2_norm();

	std::stringstream fileName;
	fileName << getenv("NATRIUM_HOME") << "/convergence-advection-solver/Level_"
			<< refinementLevel << "_p_" << fe_order << ".sp";
	std::ofstream sp_file(fileName.str().c_str());
	//streaming.getBlockSparsityPattern().print_gnuplot(sp_file);
	return result;

} /* oneTest */

} /*namespace  AdvectionBenchmark */
} /* namespace natrium */
