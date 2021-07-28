/**
 * @file SemiLagrangian_test.cpp
 * @short
 * @date 02.05.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/advection/SemiLagrangian.h"

#include "boost/test/unit_test.hpp"

/*#include "deal.II/numerics/data_out.h"
 #include "deal.II/fe/fe_dgq.h"
 #include "deal.II/fe/fe_update_flags.h"
 #include "deal.II/dofs/dof_tools.h"
 #include "deal.II/base/quadrature_lib.h"
 #include "deal.II/lac/lapack_full_matrix.h"*/

#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q19.h"

#include "natrium/advection/SemiLagrangianTools.h"

#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/PeriodicTestDomain3D.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(SemiLagrangian_test)

BOOST_AUTO_TEST_CASE(SemiLagrangian_Construction_test) {

	pout << "SemiLagrangian_Construction_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 2;
	PeriodicTestDomain2D periodic(refinementLevel);
	BOOST_CHECK_NO_THROW(
			SemiLagrangian<2> streaming(periodic, fe_order, boost::make_shared<D2Q9>(), 0.001));

	pout << "done." << endl;
} /* SemiLagrangian_Construction_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian2D_Neighborhood_test) {
	pout << "SemiLagrangian2D_Neighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic, fe_order, boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
	sl.setDeltaT(0.1);
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	// make sure cell is not in the ghost layer (begin_active() might point to a ghost cell)
	while (not cell->is_locally_owned()) {
		cell++;
	}
	Neighborhood<2> neighbors;
	sl.getNeighborhood(cell, neighbors);
	dealii::Point<2> cell_center = cell->barycenter();
	BOOST_CHECK_EQUAL(neighbors.size(), size_t(9));

	// the following tests are only valid if the cell is a left upper or left lower corner
	if (is_MPI_rank_0()) {
		double dx = pow(0.5, 3);
		size_t n_periodic_x = 0;
		size_t n_periodic_y = 0;
		for (size_t i = 0; i < neighbors.size(); i++) {
			dealii::Point<2> n_center = neighbors.at(i)->barycenter();
			// periodic boundaries
			if (n_center(0) > 0.8) {
				n_center(0) = 1.0 - n_center(0);
				n_periodic_x++;
			}
			if (n_center(1) > 0.8) {
				n_center(1) = 1.0 - n_center(1);
				n_periodic_y++;
			}
			double dist = cell_center.distance(n_center);
			// check dist \in {0, dx, sqrt(2)*dx}
			BOOST_CHECK_LE(dist, sqrt(2) * dx + 1e-8);
			if (dist < dx - 1e-8) {
				BOOST_CHECK_SMALL(dist, 1e-8);
			} else if (dist < dx + 1e-8) {
				BOOST_CHECK_GT(dist, dx - 1e-8);
			} else { // (dist < dx * sqrt(2) + 1e-8)
				BOOST_CHECK_GT(dist, dx * sqrt(2) - 1e-8);
			}
		}

		// cell is corner!
		BOOST_CHECK_EQUAL(n_periodic_x, size_t(3));
		BOOST_CHECK_EQUAL(n_periodic_y, size_t(3));
	}

	pout << "done." << endl;
} /* SemiLagrangian2D_Neighborhood_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian3D_Neighborhood_test) {
	pout << "SemiLagrangian3D_Neighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain3D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<3> sl(periodic, fe_order, boost::make_shared<D3Q19>(),
			0.001);
	sl.setupDoFs();
	sl.setDeltaT(0.1);

        typename dealii::DoFHandler<3>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	// make sure cell is not in the ghost layer
	while (not cell->is_locally_owned()) {
		cell++;
	}
	Neighborhood<3> neighbors;
	sl.getNeighborhood(cell, neighbors);
	dealii::Point<3> cell_center = cell->barycenter();
	BOOST_CHECK_EQUAL(neighbors.size(), size_t(27));

	// the following tests are only valid if the cell is a left upper or left lower corner
	if (is_MPI_rank_0()) {
		double dx = pow(0.5, 3);
		size_t n_periodic_x = 0;
		size_t n_periodic_y = 0;
		size_t n_periodic_z = 0;
		for (size_t i = 0; i < neighbors.size(); i++) {
			dealii::Point<3> n_center = neighbors.at(i)->barycenter();
			// periodic boundaries
			if (n_center(0) > 0.5) {
				n_center(0) = 1.0 - n_center(0);
				n_periodic_x++;
			}
			if (n_center(1) > 0.5) {
				n_center(1) = 1.0 - n_center(1);
				n_periodic_y++;
			}
			if (n_center(2) > 0.5) {
				n_center(2) = 1.0 - n_center(2);
				n_periodic_z++;
			}
			double dist = cell_center.distance(n_center);
			// check dist \in {0, dx, sqrt(2)*dx sqrt(3)*dx}
			BOOST_CHECK_LE(dist, sqrt(3) * dx + 1e-8);
			if (dist < dx - 1e-8) {
				BOOST_CHECK_SMALL(dist, 1e-8);
			} else if (dist < dx + 1e-8) {
				BOOST_CHECK_GT(dist, dx - 1e-8);
			} else if (dist < dx * sqrt(2) + 1e-8) {
				BOOST_CHECK_GT(dist, dx * sqrt(2) - 1e-8);
			} else { // (dist < dx * sqrt(3) + 1e-8)
				BOOST_CHECK_GT(dist, dx * sqrt(3) - 1e-8);
			}
		}

		// cell is corner!
		BOOST_CHECK_EQUAL(n_periodic_x, size_t(9));
		BOOST_CHECK_EQUAL(n_periodic_y, size_t(9));
		BOOST_CHECK_EQUAL(n_periodic_z, size_t(9));
	}

	pout << "done." << endl;
} /* SemiLagrangian3D_Neighborhood_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian2D_FindPointInNeighborhood_test) {
	pout << "SemiLagrangian2D_FindPointInNeighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic, fe_order, boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
	sl.setDeltaT(0.1);

        typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();

	dealii::Point<2> p(1, 2);
	typename dealii::DoFHandler<2>::active_cell_iterator found_in =
			sl.recursivelySearchInNeighborhood(p, cell);
	BOOST_CHECK(found_in == sl.getDoFHandler()->end());

	dealii::Point<2> p2(1, 0.5);
	found_in = sl.recursivelySearchInNeighborhood(p2, cell);
	BOOST_CHECK(found_in != sl.getDoFHandler()->end());

	pout << "done." << endl;
} /* SemiLagrangian2D_FindPointInNeighborhood_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian3D_FindPointInNeighborhood_test) {
	pout << "SemiLagrangian3D_FindPointInNeighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain3D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<3> sl(periodic, fe_order, boost::make_shared<D3Q19>(),
			0.001);
	sl.setupDoFs();
	sl.setDeltaT(0.1);

        typename dealii::DoFHandler<3>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();

	dealii::Point<3> p(1, 2, 1);
	typename dealii::DoFHandler<3>::active_cell_iterator found_in =
			sl.recursivelySearchInNeighborhood(p, cell);
	if (is_MPI_rank_0())
		BOOST_CHECK(found_in == sl.getDoFHandler()->end());

	dealii::Point<3> p2(1, 0.5, 0.7);
	found_in = sl.recursivelySearchInNeighborhood(p2, cell);
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1) {
		// otherwise, p2 could belong to another process' domain
		BOOST_CHECK(found_in != sl.getDoFHandler()->end());
	}

	pout << "done." << endl;
} /* SemiLagrangian3D_FindPointInNeighborhood_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian2D_FaceCrossedFirst_test) {
	pout << "SemiLagrangian2D_FaceCrossedFirst_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic, fe_order, boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	// make sure cell is not in the ghost layer (begin_active() might point to a ghost cell)
	while (not cell->is_locally_owned()) {
		cell++;
	}

	for (size_t without_mapping = 0; without_mapping <= 1; without_mapping++) {
		if (without_mapping)
			pout << "... without mapping" << endl;
		else
			pout << "... with mapping" << endl;
		// the following tests are only valid if the cell is the left lower corner
		if (is_MPI_rank_0()) {
			BOOST_CHECK(cell->at_boundary(0));
			BOOST_CHECK(cell->at_boundary(2));

			// p = (0.0625, 0.0625)
			dealii::Point<2> p = cell->barycenter();

			double lambda;
			size_t child_id;
			dealii::Point<2> pb;
			// no face crossed
			dealii::Point<2> p2(0.0625, 0.0625);
			BOOST_CHECK_EQUAL(
					sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
							without_mapping), -1);
			// face 0 crossed
			p2[0] = -0.0625;
			int result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda,
					&child_id, without_mapping);
			BOOST_CHECK_EQUAL(result, 0);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_SMALL(pb[0], 1e-20);
			BOOST_CHECK_CLOSE(pb[1], 0.0625, 1e-20);
			// face 2 crossed
			p2[0] = 0.0625;
			p2[1] = -0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 2);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.0625, 1e-20);
			BOOST_CHECK_SMALL(pb[1], 1e-20);
			// face 0 and 2 crossed, but 0 first
			p2[0] = 0.0625 - 0.25;
			p2[1] = -0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 0);
			BOOST_CHECK_CLOSE(lambda, 0.25, 1e-20);
			BOOST_CHECK_SMALL(pb[0], 1e-20);
			BOOST_CHECK_CLOSE(pb[1], 0.03125, 1e-20);
			// face 0 and 2 crossed, but 2 first
			p2[0] = -0.0625;
			p2[1] = 0.0625 - 0.25;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 2);
			BOOST_CHECK_CLOSE(lambda, 0.25, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.03125, 1e-20);
			BOOST_CHECK_SMALL(pb[1], 1e-20);
			// face 1 crossed
			p2[0] = 0.1875;
			p2[1] = 0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 1);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.125, 1e-20);
			BOOST_CHECK_CLOSE(pb[1], 0.0625, 1e-20);
			// face 3 crossed
			p2[0] = 0.0625;
			p2[1] = 0.1875;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 3);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.0625, 1e-20);
			BOOST_CHECK_CLOSE(pb[1], 0.125, 1e-20);
		}
	} /* with/without mapping */

	pout << "done." << endl;
} /* SemiLagrangian2D_FaceCrossedFirst_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian3D_FaceCrossedFirst_test) {
	pout << "SemiLagrangian3D_FaceCrossedFirst_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain3D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<3> sl(periodic, fe_order, boost::make_shared<D3Q19>(),
			0.001);
	sl.setupDoFs();
	typename dealii::DoFHandler<3>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	// make sure cell is not in the ghost layer (begin_active() might point to a ghost cell)
	/*while (not cell->is_locally_owned()) {
	 cell++;
	 }*/

	for (size_t without_mapping = 0; without_mapping < 1; without_mapping++) {
		if (without_mapping)
			pout << "... without mapping" << endl;
		else
			pout << "... with mapping" << endl;
		// the following tests are only valid if the cell is the left lower corner
		if (is_MPI_rank_0()) {
			BOOST_CHECK(cell->at_boundary(0));
			BOOST_CHECK(cell->at_boundary(2));
			BOOST_CHECK(cell->at_boundary(4));

			// p = (0.0625, 0.0625, 0.0625)
			dealii::Point<3> p = cell->barycenter();

			double lambda;
			size_t child_id;
			dealii::Point<3> pb;
			// no face crossed
			dealii::Point<3> p2(0.0625, 0.0625, 0.0625);
			BOOST_CHECK_EQUAL(
					sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
							without_mapping), -1);
			// face 0 crossed
			p2[0] = -0.0625;
			int result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda,
					&child_id, without_mapping);
			BOOST_CHECK_EQUAL(result, 0);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_SMALL(pb[0], 1e-20);
			BOOST_CHECK_CLOSE(pb[1], 0.0625, 1e-20);
			BOOST_CHECK_CLOSE(pb[2], 0.0625, 1e-20);
			// face 1 crossed
			p2[0] = 0.125 + 0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 1);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.125, 1e-13);
			BOOST_CHECK_CLOSE(pb[1], 0.0625, 1e-13);
			BOOST_CHECK_CLOSE(pb[2], 0.0625, 1e-13);
			// face 2 crossed
			p2[0] = 0.0625;
			p2[1] = -0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 2);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.0625, 1e-13);
			BOOST_CHECK_SMALL(pb[1], 1e-13);
			BOOST_CHECK_CLOSE(pb[2], 0.0625, 1e-13);
			// face 3 crossed
			p2[1] = 0.125 + 0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 3);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.0625, 1e-13);
			BOOST_CHECK_CLOSE(pb[1], 0.125, 1e-13);
			BOOST_CHECK_CLOSE(pb[2], 0.0625, 1e-13);
			// face 4 crossed
			p2[1] = 0.0625;
			p2[2] = -0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 4);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-13);
			BOOST_CHECK_CLOSE(pb[0], 0.0625, 1e-13);
			BOOST_CHECK_CLOSE(pb[1], 0.0625, 1e-13);
			BOOST_CHECK_SMALL(pb[2], 1e-13);
			// face 5 crossed
			p2[2] = 0.125 + 0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 5);
			BOOST_CHECK_CLOSE(lambda, 0.5, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.0625, 1e-13);
			BOOST_CHECK_CLOSE(pb[1], 0.0625, 1e-13);
			BOOST_CHECK_CLOSE(pb[2], 0.125, 1e-13);

			// face 0, 2 and 4 crossed, but 0 first
			p2[0] = 0.0625 - 0.25;
			p2[1] = -0.0625;
			p2[2] = -0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 0);
			BOOST_CHECK_CLOSE(lambda, 0.25, 1e-20);
			BOOST_CHECK_SMALL(pb[0], 1e-13);
			BOOST_CHECK_CLOSE(pb[1], 0.03125, 1e-13);
			BOOST_CHECK_CLOSE(pb[2], 0.03125, 1e-13);
			// face 0, 2 and 4 crossed, but 2 first
			p2[0] = -0.0625;
			p2[1] = 0.0625 - 0.25;
			p2[2] = -0.0625;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 2);
			BOOST_CHECK_CLOSE(lambda, 0.25, 1e-20);
			BOOST_CHECK_CLOSE(pb[0], 0.03125, 1e-13);
			BOOST_CHECK_SMALL(pb[1], 1e-13);
			BOOST_CHECK_CLOSE(pb[2], 0.03125, 1e-13);
			// face 0, 2 and 4 crossed, but 4 first
			p2[0] = -0.0625;
			p2[1] = -0.0625;
			p2[2] = 0.0625 - 0.25;
			result = sl.faceCrossedFirst(cell, p, p2, pb, &lambda, &child_id,
					without_mapping);
			BOOST_CHECK_EQUAL(result, 4);
			BOOST_CHECK_CLOSE(lambda, 0.25, 1e-13);
			BOOST_CHECK_CLOSE(pb[0], 0.03125, 1e-13);
			BOOST_CHECK_CLOSE(pb[1], 0.03125, 1e-13);
			BOOST_CHECK_SMALL(pb[2], 1e-13);
		}
	} /* with/without mapping */

	pout << "done." << endl;
} /* SemiLagrangian3D_FaceCrossedFirst_test */

//BOOST_AUTO_TEST_CASE(SemiLagrangian2D_SparsityPattern_test) {
//	pout << "SemiLagrangian2D_SparsityPattern_test..." << endl;
//
//	// setup system
//	size_t fe_order = 2;
//	size_t refinementLevel = 3;
//
//	PeriodicTestDomain2D periodic(refinementLevel);
//	periodic.refineAndTransform();
//	SemiLagrangian<2> sl(periodic, fe_order, boost::make_shared<D2Q9>(), 0.001);
//	sl.setupDoFs();
//
//	// check number of corresponding dofs for each dof
//	// it has to be = (fe_order + 1)^dim
//	const std::vector<std::vector<dealii::TrilinosWrappers::SparsityPattern> > & sp =
//			sl.getBlockSparsityPattern();
//	BOOST_CHECK_EQUAL(sp.size(), size_t(8));
//	for (size_t i = 0; i < sp.size(); i++) {
//		BOOST_CHECK_EQUAL(sp[i].size(), size_t(8));
//		/* the following test worked for the suboptimal sparsity patterns
//		 BOOST_CHECK_SMALL(
//		 sp[i][i].n_nonzero_elements()
//		 - sl.getDoFHandler()->n_dofs() * pow((fe_order + 1), 2),
//		 2.0);
//		 */
//	}
//
//	pout << "done." << endl;
//} /* SemiLagrangian2D_SparsityPattern_test */
//
//BOOST_AUTO_TEST_CASE(SemiLagrangian3D_SparsityPattern_test) {
//	pout << "SemiLagrangian3D_SparsityPattern_test..." << endl;
//
//	// setup system
//	size_t fe_order = 1;
//	size_t refinementLevel = 3;
//
//	PeriodicTestDomain3D periodic(refinementLevel);
//	periodic.refineAndTransform();
//	SemiLagrangian<3> sl(periodic, fe_order, boost::make_shared<D3Q19>(),
//			0.001);
//	sl.setupDoFs();
//
//	// check number of corresponding dofs for each dof
//	// it has to be = (fe_order + 1)^dim
//	const std::vector<std::vector<dealii::TrilinosWrappers::SparsityPattern> > & sp =
//			sl.getBlockSparsityPattern();
//	BOOST_CHECK_EQUAL(sp.size(), size_t(18));
//	for (size_t i = 0; i < sp.size(); i++) {
//		BOOST_CHECK_EQUAL(sp[i].size(), size_t(18));
//		/* the following test worked for the suboptimal sparsity patterns
//		 BOOST_CHECK_SMALL(
//		 sp[i][i].n_nonzero_elements()
//		 - sl.getDoFHandler()->n_dofs() * pow((fe_order + 1), 3),
//		 2.0);
//		 */
//	}
//
//	pout << "done." << endl;
//} /* SemiLagrangian3D_SparsityPattern_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian2D_ConstantStreaming_test) {
	pout << "SemiLagrangian2D_ConstantStreaming_test..." << endl;

	// setup system
	size_t fe_order = 3;
	size_t refinementLevel = 3;

	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic, fe_order, boost::make_shared<D2Q9>(), 0.01);
	sl.setupDoFs();
    sl.setDeltaT(0.1);

    sl.reassemble();

	distributed_block_vector ones;
	distributed_block_vector result;
	ones.reinit(8);
	result.reinit(8);
	for (size_t i = 0; i < 8; i++) {
		ones.block(i).reinit(sl.getLocallyOwnedDofs(), MPI_COMM_WORLD);
		result.block(i).reinit(sl.getLocallyOwnedDofs(), MPI_COMM_WORLD);
		// reinit does only change the size but not the content
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(sl.getLocallyOwnedDofs().begin());
		dealii::IndexSet::ElementIterator end(sl.getLocallyOwnedDofs().end());
		for (; it != end; it++) {
			size_t j = *it;
			ones.block(i)(j) = 1;
		}
	}

	//sl.getSystemMatrix().print(cout);
	sl.getSystemMatrix().vmult(result, ones);
	result -= ones;
	BOOST_CHECK_LE(result.norm_sqr(), 1e-6);

	pout << "done." << endl;
} /* SemiLagrangian2D_ConstantStreaming_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian3D_ConstantStreaming_test) {
	pout << "SemiLagrangian3D_ConstantStreaming_test..." << endl;

	// setup system
	size_t fe_order = 1;
	size_t refinementLevel = 3;

	PeriodicTestDomain3D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<3> sl(periodic, fe_order, boost::make_shared<D3Q19>(), 0.01);
	sl.setupDoFs();
    sl.setDeltaT(0.1);

    sl.reassemble();

	distributed_block_vector ones;
	distributed_block_vector result;
	ones.reinit(18);
	result.reinit(18);
	for (size_t i = 0; i < 18; i++) {
		ones.block(i).reinit(sl.getLocallyOwnedDofs(), MPI_COMM_WORLD);
		result.block(i).reinit(sl.getLocallyOwnedDofs(), MPI_COMM_WORLD);
		// reinit does only change the size but not the content
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(sl.getLocallyOwnedDofs().begin());
		dealii::IndexSet::ElementIterator end(sl.getLocallyOwnedDofs().end());
		for (; it != end; it++) {
			size_t j = *it;
			ones.block(i)(j) = 1;
		}
	}

	sl.getSystemMatrix().vmult(result, ones);
	result -= ones;
	BOOST_CHECK_LE(result.norm_sqr(), 1e-6);

	pout << "done." << endl;
} /* SemiLagrangian3D_ConstantStreaming_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian3D_VectorReference_test) {
	pout << "SemiLagrangian3D_VectorReference_test..." << endl;

	// setup system
	size_t fe_order = 1;
	size_t refinementLevel = 3;

	PeriodicTestDomain3D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<3> sl(periodic, fe_order, boost::make_shared<D3Q19>(), 0.01);
	sl.setupDoFs();
    sl.setDeltaT(0.1);

    sl.reassemble();

	distributed_block_vector ones;
	distributed_block_vector result;
	ones.reinit(18);
	result.reinit(18);
	for (size_t i = 0; i < 18; i++) {
		ones.block(i).reinit(sl.getLocallyOwnedDofs(), MPI_COMM_WORLD);
		result.block(i).reinit(sl.getLocallyOwnedDofs(), MPI_COMM_WORLD);
		// reinit does only change the size but not the content
		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(sl.getLocallyOwnedDofs().begin());
		dealii::IndexSet::ElementIterator end(sl.getLocallyOwnedDofs().end());
		for (; it != end; it++) {
			size_t j = *it;
			ones.block(i)(j) = 1;
		}
	}

	pout << "done." << endl;
} /* SemiLagrangian3D_VectorReference_test */

BOOST_AUTO_TEST_SUITE_END()
