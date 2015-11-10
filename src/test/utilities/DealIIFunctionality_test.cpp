/**
 * @file DealIIFunctionality_test.cpp
 * @short
 * @date 07.10.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/BasicNames.h"

#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_generator.h"

#include "boost/test/unit_test.hpp"

namespace natrium {

BOOST_AUTO_TEST_SUITE(DealIIFunctionality_test)

BOOST_AUTO_TEST_CASE(DealIIFunctionality_Periodicity_test) {
	// Check how the built-in periodic bc works in deal.II
	pout << "DealIIFunctionality_Periodicity_test..." << endl;

	Mesh<2> square(MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_cube(square, 0, 1);
	Mesh<2>::active_cell_iterator cell = square.begin_active();
	Mesh<2>::active_cell_iterator neighbor = square.begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom


	// check if periodicity updates neighbors
	BOOST_CHECK_EQUAL(cell->neighbor_index(0), -1);
	std::vector<
			dealii::GridTools::PeriodicFacePair<typename Mesh<2>::cell_iterator> > matched_pairs;
	dealii::GridTools::collect_periodic_faces(square, 0, 1, 0, matched_pairs);
	square.add_periodicity(matched_pairs);
	// BOOST_CHECK_NE(cell->neighbor_index(0), -1); ->FALSE
	// So: add_periodicity does not update neighbor relations, just the ghost layer

	// check if cell is still at boundary
	BOOST_CHECK(cell->at_boundary());

	// check if periodicity is inherited by subcells
	square.refine_global(1);

	// check if neighbor relations work generally
	// therefore pick upper right corner cell
	for (cell = square.begin_active(); cell != square.end(); cell ++)
		if (cell->at_boundary(2))
			if (cell->at_boundary(1))
				break;
	BOOST_CHECK_NE(cell->neighbor_index(0), -1);

	// check if neighbor relations can be enforced through set_neighbor
	//----------------------------------------------------------------
	// ---------------------------------------------------------------
	//
	//      |                   |                       |
	//      |     neighbor      |           cell        |
	//      |                   |                       |
	//
	//
	// neighbor->face(0)            				cell->face(1)
	//----------------------------------------------------------------
	//----------------------------------------------------------------

	for (cell = square.begin_active(); cell != square.end(); cell ++)
		if (cell->at_boundary(2))
			if (cell->at_boundary(1))
				break;
	for (neighbor = square.begin_active(); neighbor != square.end(); neighbor ++)
		if (neighbor->at_boundary(2))
			if (neighbor->at_boundary(0))
				break;
	cell->set_neighbor(1, neighbor);
	neighbor->set_neighbor(0, cell);
	BOOST_CHECK_NO_THROW(cell->neighbor(1));
	BOOST_CHECK(cell->neighbor(1) == neighbor);
	BOOST_CHECK(neighbor->neighbor(0) == cell);
	//BOOST_CHECK(cell->neighbor_is_coarser(1) == false);
	// this fails: only  cell->neighbor(i)  gives the right thing , but the more sophisticated functions do not work
	// consequently, set_neighbor is not sufficient
	//BOOST_CHECK(cell->neighbor_of_neighbor(1) == 0);
	//BOOST_CHECK(neighbor->neighbor_of_neighbor(0) == 1);




	pout << "done." << endl;

}


BOOST_AUTO_TEST_CASE(DealIIFunctionality_UserFlags_test) {
	// Check how the built-in periodic bc works in deal.II
	pout << "DealIIFunctionality_UserFlags_test..." << endl;

	Mesh<2> square(MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_cube(square, 0, 1);
	Mesh<2>::active_cell_iterator cell = square.begin_active();
	cell->face(0)->set_user_flag();

	BOOST_CHECK(cell->face(0)->user_flag_set());

	// setting twice -> still true?
	cell->face(0)->set_user_flag();
	BOOST_CHECK(cell->face(0)->user_flag_set());

	// different iterator -> still true?
	Mesh<2>::cell_iterator cell2 = square.begin();
	BOOST_CHECK(cell2->face(0)->user_flag_set());

	pout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
