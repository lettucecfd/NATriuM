/*
 * Checkpoint_test.cpp
 *
 *  Created on: 19.04.2016
 *      Author: akraem3m
 */

#include "natrium/solver/Checkpoint.h"

#include "boost/test/unit_test.hpp"
#include "boost/filesystem.hpp"

#include "deal.II/base/index_set.h"
#include "deal.II/base/function_lib.h"
#include "deal.II/numerics/vector_tools.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/benchmarks/TaylorGreenVortex2D.h"
#include "natrium/advection/SEDGMinLee.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/solver/CFDSolver.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(Checkpoint_test)

BOOST_AUTO_TEST_CASE(Checkpoint_Create_test) {
	pout << "Checkpoint_Create_test..." << endl;

	// prepare directory
	boost::filesystem::path check_dir = "/tmp/natrium_checkpoint_test";
	if (is_MPI_rank_0()) {
		boost::filesystem::create_directory(check_dir);
	}

	//test
	Checkpoint<2> c(0, check_dir);

	// clean up
	if (is_MPI_rank_0()) {
		boost::filesystem::remove_all(check_dir);
	}
	pout << "done" << endl;
} /* Checkpoint_Create_test */

BOOST_AUTO_TEST_CASE(Checkpoint_SaveAndLoadUnchanged_test) {
	pout << "Checkpoint_SaveAndLoadUnchanged_test..." << endl;

	// for sedg and sl streaming
	for (size_t sl = 0; sl <= 1; sl++) {
		if (sl) {
			pout << " - for semi-Lagrangian streaming" << endl;
		} else {
			pout << " - for SEDG streaming" << endl;
		}

		// prepare directory
		boost::filesystem::path check_dir = "/tmp/natrium_checkpoint_test";
		if (is_MPI_rank_0()) {
			boost::filesystem::create_directory(check_dir);
		}

		// prepare problem and advection operator
		size_t ref_level = 3;
		double viscosity = 1;
		size_t p = 2;
		size_t Q = 9;
		boost::shared_ptr<ProblemDescription<2> > tgv = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, ref_level);
		tgv->refineAndTransform();

		boost::shared_ptr<AdvectionOperator<2> > advection;
		if (sl) {
			advection = boost::make_shared<SemiLagrangian<2> >(*tgv, p, QGAUSS,
					GAUSS_LOBATTO_CHEBYSHEV_POINTS, boost::make_shared<D2Q9>(1),
					0.1);
		} else {
			advection = boost::make_shared<SEDGMinLee<2> >(*tgv, p,
					boost::make_shared<D2Q9>(1));
		}
		advection->setupDoFs();
		DistributionFunctions f(not sl);
		//if (sl)
			f.reinit(Q, advection->getLocallyOwnedDofs(),
					advection->getLocallyRelevantDofs(), MPI_COMM_WORLD, not sl);

		//else
		//	f.reinit(Q, advection->getLocallyOwnedDofs(), MPI_COMM_WORLD);

		// random initialization
		for (size_t i = 0; i < Q; i++) {
			dealii::IndexSet::ElementIterator it(
					advection->getLocallyOwnedDofs().begin());
			dealii::IndexSet::ElementIterator end(
					advection->getLocallyOwnedDofs().end());
			for (it = advection->getLocallyOwnedDofs().begin(); it != end;
					it++) {
				size_t j = *it;
				f.at(i)(j) = random();
			}
		}
		f.compress(dealii::VectorOperation::insert);

		// write checkpoint
		CheckpointStatus status;
		status.feOrder = p;
		status.iterationNumber = 1;
		status.stencilScaling = 1;
		status.time = 10.0;
		Checkpoint<2> c1(1, check_dir);
		c1.write(*tgv->getMesh(), f, *advection->getDoFHandler(), status);

		// read checkpoint
		DistributionFunctions f2(not sl);
		CheckpointStatus status2;
		boost::shared_ptr<ProblemDescription<2> > tgv2 = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, ref_level);
		boost::shared_ptr<AdvectionOperator<2> > advection2;
		if (sl) {
			advection2 = boost::make_shared<SemiLagrangian<2> >(*tgv2, p,
					QGAUSS, GAUSS_LOBATTO_CHEBYSHEV_POINTS,
					boost::make_shared<D2Q9>(1), 0.1);
		} else {
			advection2 = boost::make_shared<SEDGMinLee<2> >(*tgv2, p,
					boost::make_shared<D2Q9>(1));
		}
		Checkpoint<2> c2(1, check_dir);
		c2.load(f2, *tgv2, *advection2, status2);

		// compare
		BOOST_CHECK_EQUAL(status2.feOrder, size_t(p));
		BOOST_CHECK_EQUAL(status2.iterationNumber, size_t(1));
		BOOST_CHECK_CLOSE(status2.stencilScaling, 1.0, 1e-8);
		BOOST_CHECK_CLOSE(status2.time, 10.0, 1e-8);
		BOOST_CHECK_EQUAL(f2.size(), Q);
		//BOOST_CHECK(f2.equals(f));
		DistributionFunctions diff(f2);
		for (size_t i = 0; i < Q; i++) {
			diff.at(i) *= -1.0;
			diff.at(i).add(f.at(i), true);
			BOOST_CHECK_SMALL(diff.at(i).norm_sqr(), 1e-8);
		}

		// clean up
		if (is_MPI_rank_0()) {
			boost::filesystem::remove_all(check_dir);
		}
	} /* end for sl and sedg */

	pout << "done" << endl;
} /* Checkpoint_SaveAndLoadUnchanged_test */

BOOST_AUTO_TEST_CASE(Checkpoint_ResumeRefined_test) {
	pout << "Checkpoint_ResumeRefined_test..." << endl;

	// for sedg and sl streaming
	for (size_t sl = 0; sl <= 1; sl++) {
		if (sl) {
			pout << " - for semi-Lagrangian streaming" << endl;
		} else {
			pout << " - for SEDG streaming" << endl;
		}

// prepare directory
		boost::filesystem::path check_dir = "/tmp/natrium_checkpoint_test";
		if (is_MPI_rank_0()) {
			boost::filesystem::create_directory(check_dir);
		}

// prepare problem and advection operator
		size_t ref_level = 3;
		double viscosity = 1;
		size_t p = 2;
		size_t Q = 9;
		boost::shared_ptr<ProblemDescription<2> > tgv = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, ref_level);
		tgv->refineAndTransform();

		boost::shared_ptr<AdvectionOperator<2> > advection;
		if (sl) {
			advection = boost::make_shared<SemiLagrangian<2> >(*tgv, p, QGAUSS,
					GAUSS_LOBATTO_CHEBYSHEV_POINTS, boost::make_shared<D2Q9>(1),
					0.1);
		} else {
			advection = boost::make_shared<SEDGMinLee<2> >(*tgv, p,
					boost::make_shared<D2Q9>(1));
		}
		advection->setupDoFs();
		DistributionFunctions f(not sl);
		//if (sl)
			f.reinit(Q, advection->getLocallyOwnedDofs(),
					advection->getLocallyRelevantDofs(), MPI_COMM_WORLD, not sl);
		//else
		//	f.reinit(Q, advection->getLocallyOwnedDofs(), MPI_COMM_WORLD);

// initialization
		dealii::Point<2> exponents(2, 1);
		dealii::Functions::Monomial<2> mon(exponents); // has to be interpolated exactly
		for (size_t i = 0; i < Q; i++) {
			dealii::VectorTools::interpolate(advection->getMapping(),
					*advection->getDoFHandler(), mon, f.at(i));
		}

// write checkpoint
		CheckpointStatus status;
		status.feOrder = p;
		status.iterationNumber = 1;
		status.stencilScaling = 1;
		status.time = 10.0;
		Checkpoint<2> c1(1, check_dir);
		c1.write(*tgv->getMesh(), f, *advection->getDoFHandler(), status);

// read checkpoint with different ref_level
		ref_level = 4;
		DistributionFunctions f2(not sl);
		CheckpointStatus status2;
		boost::shared_ptr<ProblemDescription<2> > tgv2 = boost::make_shared<
				TaylorGreenVortex2D>(viscosity, ref_level);
		boost::shared_ptr<AdvectionOperator<2> > advection2;
		if (sl) {
			advection2 = boost::make_shared<SemiLagrangian<2> >(*tgv2, p,
					QGAUSS, GAUSS_LOBATTO_CHEBYSHEV_POINTS,
					boost::make_shared<D2Q9>(1), 0.1);
		} else {
			advection2 = boost::make_shared<SEDGMinLee<2> >(*tgv2, p,
					boost::make_shared<D2Q9>(1));
		}
		Checkpoint<2> c2(1, check_dir);
		c2.load(f2, *tgv2, *advection2, status2);

// compare
		BOOST_CHECK_EQUAL(status2.feOrder, size_t(p));
		BOOST_CHECK_EQUAL(status2.iterationNumber, size_t(1));
		BOOST_CHECK_CLOSE(status2.stencilScaling, 1.0, 1e-8);
		BOOST_CHECK_CLOSE(status2.time, 10.0, 1e-8);
		BOOST_CHECK_EQUAL(f2.size(), Q);
// check interpolation
		distributed_vector expected;
		expected.reinit(advection2->getLocallyOwnedDofs(), MPI_COMM_WORLD);
		dealii::VectorTools::interpolate(advection2->getMapping(),
				*advection2->getDoFHandler(), mon, expected);
		for (size_t i = 0; i < Q; i++) {
			distributed_vector diff = expected;
			diff *= -1.0;
			diff.add(f2.at(i), true); // allow different maps
			BOOST_CHECK_SMALL(diff.norm_sqr(), 1e-6);
		}

// clean up
		if (is_MPI_rank_0()) {
			boost::filesystem::remove_all(check_dir);
		}

	} /* end for sl and sedg */
	pout << "done" << endl;
} /* Checkpoint_ResumeRefined_test */

BOOST_AUTO_TEST_CASE(Checkpoint_ResumeOtherP_test) {
// Not supported by the checkpoint routine so far
} /* Checkpoint_ResumeOtherP_test */

BOOST_AUTO_TEST_CASE(Checkpoint_ResumeOtherMa) {
	pout << "Checkpoint_ResumeOtherMa..." << endl;

// prepare directory
	boost::filesystem::path check_dir = "/tmp/natrium_checkpoint_test";
	if (is_MPI_rank_0()) {
		boost::filesystem::create_directory(check_dir);
	}

// prepare and write checkpoint
	size_t ref_level = 3;
	double viscosity = 1;
	boost::shared_ptr<ProblemDescription<2> > tgv = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, ref_level);
	boost::shared_ptr<SolverConfiguration> cfg = boost::make_shared<
			SolverConfiguration>();
	cfg->setCommandLineVerbosity(ERROR);
	cfg->setUserInteraction(false);
	cfg->setOutputCheckpointInterval(1);
	cfg->setNumberOfTimeSteps(1);
	CFDSolver<2> solver(cfg, tgv);
	solver.run();
	const vector<distributed_vector>& u1 = solver.getVelocity();
	const distributed_vector& rho1 = solver.getDensity();

// load checkpoint
	boost::shared_ptr<ProblemDescription<2> > tgv2 = boost::make_shared<
			TaylorGreenVortex2D>(viscosity, ref_level);
	cfg->setStencilScaling(110.5);
	cfg->setRestartAtIteration(1);
	CFDSolver<2> solver2(cfg, tgv2);
	solver2.collide();
	const vector<distributed_vector>& u2 = solver2.getVelocity();
	const distributed_vector& rho2 = solver2.getDensity();

// compare
	distributed_vector cmp;
	CFDSolverUtilities::getWriteableDensity(cmp, rho1,
			solver2.getAdvectionOperator()->getLocallyOwnedDofs());
	distributed_vector rho2_writeable;
	vector<distributed_vector> u2_writeable;
	CFDSolverUtilities::getWriteableDensity(rho2_writeable, rho2,
			solver2.getAdvectionOperator()->getLocallyOwnedDofs());
	CFDSolverUtilities::getWriteableVelocity(u2_writeable, u2,
			solver2.getAdvectionOperator()->getLocallyOwnedDofs());
	cmp.add(-1.0, rho2_writeable);
	BOOST_CHECK_SMALL(cmp.norm_sqr(), 1e-4);
	cmp = u1.at(0);
	cmp.add(-1.0, u2_writeable.at(0));
	BOOST_CHECK_SMALL(cmp.norm_sqr(), 1e-4);
	cmp = u1.at(1);
	cmp.add(-1.0, u2_writeable.at(1));
	BOOST_CHECK_SMALL(cmp.norm_sqr(), 1e-4);

// clean up
	if (is_MPI_rank_0()) {
		boost::filesystem::remove_all(check_dir);
	}

	pout << "done" << endl;
} /* Checkpoint_ResumeOtherMa */

BOOST_AUTO_TEST_SUITE_END()
