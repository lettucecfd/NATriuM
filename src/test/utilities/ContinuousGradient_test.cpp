/**
 * @file DealIIFunctionality_test.cpp
 * @short
 * @date 07.10.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "natrium/utilities/ContinuousBoundaryGradient.h"

#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/fe/mapping_q1.h"
#include "deal.II/lac/affine_constraints.h"

#include "boost/test/unit_test.hpp"

#include "natrium/utilities/BasicNames.h"


using namespace natrium;

// Non-differentialble function
class NonDiffFunction2D: public dealii::Function<2> {
public:
	NonDiffFunction2D(){

	}
	virtual ~NonDiffFunction2D(){

	}
	virtual double 	value (const dealii::Point< 2 > &points, const unsigned int component) const{
		assert (component == 0);
		double x = points[0];
		return abs(sin(x/(2.0*M_PI)));
	}
};

class NonDiffFunction3D: public dealii::Function<3> {
public:
	NonDiffFunction3D(){

	}
	virtual ~NonDiffFunction3D(){

	}
	virtual double 	value (const dealii::Point< 3 > &points, const unsigned int component) const{
		assert (component == 0);
		double x = points[0];
		return abs(sin(x/(2.0*M_PI)));
	}
};

class DiffFunction2D: public dealii::Function<2> {
public:
	DiffFunction2D(){

	}
	virtual ~DiffFunction2D(){

	}
	virtual double 	value (const dealii::Point< 2 > &points, const unsigned int component) const{
		assert (component == 0);
		double x = points[0];
		return sin(x/(2.0*M_PI));
	}
};

class DiffFunction3D: public dealii::Function<3> {
public:
	DiffFunction3D(){

	}
	virtual ~DiffFunction3D(){

	}
	virtual double 	value (const dealii::Point< 3 > &points, const unsigned int component) const{
		assert (component == 0);
		double x = points[0];
		return sin(x/(2.0*M_PI));
	}
};


BOOST_AUTO_TEST_SUITE(ContinuousBoundaryGradient_test)

BOOST_AUTO_TEST_CASE(ContinuousBoundaryGradient_Initialization_test) {
	pout << "ContinuousBoundaryGradient_Initialization_test..." << endl;

	// prepare mesh
	Mesh<2> square(MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_cube(square, 0, 1);
	Mesh<2>::active_cell_iterator cell = square.begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom
	square.refine_global(4);
	std::set<size_t> boundary_ids;
	boundary_ids.insert(0);
	boundary_ids.insert(2);

	// make dofs and solution
	dealii::QGaussLobatto<1> q1(4);
	dealii::QGaussLobatto<2> q2(4);
	dealii::DoFHandler<2> dof(square);
	dealii::FE_Q<2> fe (q1);
	dof.distribute_dofs(fe);
	distributed_vector u(dof.locally_owned_dofs());
	dealii::MappingQ1<2> m;
	dealii::AffineConstraints<double> c;
	c.close();
	// project does not work in parallel
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1) {
		dealii::VectorTools::project(m, dof, c, q2, DiffFunction2D(), u);
		ContinuousBoundaryGradient<2> grad(dof, m,  q1, q1, q2 , boundary_ids);
		grad.reinit();
	}

	dof.clear();
	pout << "done." << endl;

} /* ContinuousBoundaryGradient_Initialization_test */

BOOST_AUTO_TEST_CASE(ContinuousBoundaryGradient_Functionality2D_test) {
	pout << "ContinuousBoundaryGradient_Functionality2D_test..." << endl;

	// prepare mesh
	Mesh<2> square(MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_cube(square, 0, 1);
	Mesh<2>::active_cell_iterator cell0 = square.begin_active();
	cell0->face(0)->set_all_boundary_ids(0);  // left
	cell0->face(1)->set_all_boundary_ids(1);  // right
	cell0->face(2)->set_all_boundary_ids(2);  // top
	cell0->face(3)->set_all_boundary_ids(3);  // bottom
	square.refine_global(4);
	std::set<size_t> boundary_ids;
	boundary_ids.insert(2);

	// make dofs and solution
	dealii::QGaussLobatto<1> q1(4);
	dealii::QGaussLobatto<2> q2(4);
	dealii::DoFHandler<2> dof(square);
	dealii::FE_Q<2> fe (q1);
	dof.distribute_dofs(fe);
	distributed_vector u(dof.locally_owned_dofs());
	dealii::MappingQ1<2> m;
	dealii::AffineConstraints<double> c;
	c.close();
	// project does not work in parallel
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1) {
		dealii::VectorTools::project(m, dof, c, q2, DiffFunction2D(), u);
		ContinuousBoundaryGradient<2> grad(dof, m,  q1, q1, q2 , boundary_ids);
		grad.reinit();
	}


/*
	cout << "Hallo " << endl;
	grad.calculateGradients(u);

	dealii::UpdateFlags flags = dealii::update_gradients;
	dealii::FEFaceValues<2> face_values(m, dof.get_fe(),
			q1, flags);
	std::vector<dealii::types::global_dof_index> indices;
	indices.resize(dof.get_fe().dofs_per_cell);

	std::vector<dealii::Tensor<1, 2> > gradients;
	size_t n_q_points = q1.size();
	gradients.resize(n_q_points);
	typename dealii::DoFHandler<2>::active_cell_iterator cell = dof.begin_active();
	typename ContinuousBoundaryGradient<2>::active_cell_iterator grad_cell = grad.begin_active();
	typename ContinuousBoundaryGradient<2>::active_cell_iterator e = grad.end();

	// check that the gradient was not changed when starting with a continuous gradient
	for (; grad_cell != e; ++cell, ++grad_cell){
		if (cell->at_boundary()){
			for (size_t f = 0; f < dealii::GeometryInfo<2>::faces_per_cell; f++){
				if (cell->face(f)->boundary_id() == 2){
					face_values.reinit(cell, f);
					face_values.get_function_gradients(u, gradients);
					for (size_t q = 0; q < n_q_points; q++){
						for (size_t i = 0; i < 2; i++){
							BOOST_CHECK_CLOSE(gradients.at(q)[i], grad.get_gradient_component(q,i), 1e-10);
						}
					}
				}
			}
		}
	}

*/

	dof.clear();

	pout << "done." << endl;
} /* ContinuousBoundaryGradient_Functionality2D_test */

BOOST_AUTO_TEST_CASE(ContinuousBoundaryGradient_Functionality3D_test) {
	//pout << "ContinuousBoundaryGradient_Functionality3D_test..." << endl;


	//pout << "done." << endl;
} /* ContinuousBoundaryGradient_Functionality3D_test */


BOOST_AUTO_TEST_SUITE_END()
