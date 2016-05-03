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
			SemiLagrangian<2> streaming(periodic.getMesh(), periodic.getBoundaries(), fe_order, boost::make_shared<D2Q9>()));

	pout << "done." << endl;
} /* SemiLagrangian_Construction_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian2D_Neighborhood_test) {
	pout << "SemiLagrangian2D_Neighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<2> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D2Q9>());
	sl.setupDoFs();
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	SemiLagrangian<2>::Neighborhood neighbors;
	sl.getNeighborhood(cell, neighbors);
	dealii::Point<2> cell_center = cell->barycenter();
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1){
		BOOST_CHECK_EQUAL(neighbors.size(), 13);
	} else {
		BOOST_CHECK_GE(neighbors.size(), 9);
	}
	double dx = pow(0.5, 3);
	size_t n_periodic_x = 0;
	size_t n_periodic_y = 0;
	for (size_t i = 0; i < neighbors.size(); i++) {
		dealii::Point<2> n_center = neighbors.at(i)->barycenter();
		// periodic boundaries
		if (n_center(0) > 0.5) {
			n_center(0) = 1.0 - n_center(0);
			n_periodic_x++;
		}
		if (n_center(1) > 0.5) {
			n_center(1) = 1.0 - n_center(1);
			n_periodic_y++;
		}
		double dist = cell_center.distance(n_center);
		// check dist \in {0, dx, sqrt(2)*dx, 2*dx}
		BOOST_CHECK_LE(dist, 2 * dx + 1e-8);
		if (dist < dx - 1e-8) {
			BOOST_CHECK_SMALL(dist, 1e-8);
		} else if (dist < dx + 1e-8) {
			BOOST_CHECK_GT(dist, dx - 1e-8);
		} else if (dist < dx * sqrt(2) + 1e-8) {
			BOOST_CHECK_GT(dist, dx * sqrt(2) - 1e-8);
		} else {
			BOOST_CHECK_GT(dist, dx * 2 - 1e-8);
		}
	}
	// cell is corner!
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1){
		BOOST_CHECK_EQUAL(n_periodic_x, 4);
		BOOST_CHECK_EQUAL(n_periodic_y, 4);
	} else {
		BOOST_CHECK_GE(n_periodic_x, 3);
		BOOST_CHECK_GE(n_periodic_y, 3);
	}

	pout << "done." << endl;
} /* SemiLagrangian2D_Neighborhod_test */

BOOST_AUTO_TEST_CASE(SemiLagrangian3D_Neighborhood_test) {
	pout << "SemiLagrangian3D_Neighborhood_test..." << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain3D periodic(refinementLevel);
	periodic.refineAndTransform();
	SemiLagrangian<3> sl(periodic.getMesh(), periodic.getBoundaries(), fe_order,
			boost::make_shared<D3Q19>());
	sl.setupDoFs();
	typename dealii::DoFHandler<3>::active_cell_iterator cell =
			sl.getDoFHandler()->begin_active();
	SemiLagrangian<3>::Neighborhood neighbors;
	sl.getNeighborhood(cell, neighbors);
	dealii::Point<3> cell_center = cell->barycenter();
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1){
		BOOST_CHECK_EQUAL(neighbors.size(), 63);
	} else {
		BOOST_CHECK_GE(neighbors.size(), 27);
	}
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
		// check dist \in {0, dx, sqrt(2)*dx sqrt(3)*dx, 2*dx}
		BOOST_CHECK_LE(dist, 3 * dx + 1e-8);
		if (dist < dx - 1e-8) {
			BOOST_CHECK_SMALL(dist, 1e-8);
		} else if (dist < dx + 1e-8) {
			BOOST_CHECK_GT(dist, dx - 1e-8);
		} else if (dist < dx * sqrt(2) + 1e-8) {
			BOOST_CHECK_GT(dist, dx * sqrt(2) - 1e-8);
		} else if (dist < dx * sqrt(3) + 1e-8) {
			BOOST_CHECK_GT(dist, dx * sqrt(3) - 1e-8);
		} else if (dist < dx * 2 + 1e-8) {
			BOOST_CHECK_GT(dist, dx * 2 - 1e-8);
		} else if (dist < dx * sqrt(5) + 1e-8) {
			BOOST_CHECK_GT(dist, dx * sqrt(5) - 1e-8);
		} else {
			BOOST_CHECK_GT(dist, dx * 3 - 1e-8);
		}
	}
	// cell is corner!
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1){
	BOOST_CHECK_EQUAL(n_periodic_x, 19);
	BOOST_CHECK_EQUAL(n_periodic_y, 19);
	BOOST_CHECK_EQUAL(n_periodic_z, 19);
	} else {
		BOOST_CHECK_GE(n_periodic_x, 13);
		BOOST_CHECK_GE(n_periodic_y, 13);
	}


	pout << "done." << endl;
} /* SemiLagrangian3D_Neighborhod_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
