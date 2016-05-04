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

#include "natrium/utilities/BasicNames.h"

#include "natrium/benchmarks/PeriodicTestDomain3D.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(SemiLagrangian_test)

BOOST_AUTO_TEST_CASE(SemiLagrangian_Construction_test) {

	pout << "SemiLagrangian_Construction_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 2;
	PeriodicTestDomain2D periodic(refinementLevel);
	BOOST_CHECK_NO_THROW(
			SemiLagrangian<2> streaming(periodic.getMesh(), periodic.getBoundaries(), fe_order, boost::make_shared<D2Q9>(), 0.001));

	pout << "done." << endl;
} /* SemiLagrangian_Construction_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian2D_Neighborhood_test) {
	pout << "SemiLagrangian2D_Neighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	// make sure cell is not in the ghost layer (begin_active() might point to a ghost cell)
	while (not cell->is_locally_owned()) {
		cell++;
	}
	SemiLagrangian<2>::Neighborhood neighbors;
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
	SemiLagrangian<3> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D3Q19>(), 0.001);
	sl.setupDoFs();
	typename dealii::DoFHandler<3>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	// make sure cell is not in the ghost layer
	while (not cell->is_locally_owned()) {
		cell++;
	}
	SemiLagrangian<3>::Neighborhood neighbors;
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
	SemiLagrangian<2> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D2Q9>(), 0.001);
	sl.setupDoFs();
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
	SemiLagrangian<3> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D3Q19>(), 0.001);
	sl.setupDoFs();
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

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
