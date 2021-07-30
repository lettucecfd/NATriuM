/**
 * @file Cylinder2D.cpp
 * @short 
 * @date 09.10.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Cylinder2D.h"

#include <fstream>

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_in.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/DoNothingBoundary.h"
#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/DealiiExtensions.h"
#include "natrium/utilities/Logging.h"

namespace natrium {

Cylinder2D::Cylinder2D(double viscosity, double inletVelocity) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_inletVelocity(
				inletVelocity) {
	setCharacteristicLength(1.0);

	this->setInitialU(boost::make_shared<InitialU>(this));
	/// apply boundary values
	setBoundaries(makeBoundaries(inletVelocity));
}

Cylinder2D::~Cylinder2D() {
}

void Cylinder2D::make_inner_manifold(dealii::Triangulation<2>& mesh,
		const dealii::SphericalManifold<2>& manifold, double radius, int boundary_id,
		double tol) {

	dealii::Triangulation<2>::active_cell_iterator cell = mesh.begin_active(),
			endc = mesh.end();
	for (; cell != endc; ++cell) {
		for (unsigned int f = 0; f < dealii::GeometryInfo<2>::faces_per_cell;
				++f) {
			bool is_inner_face = true;
			for (unsigned int v = 0;
					v < dealii::GeometryInfo<2>::vertices_per_face; ++v) {
				const double distance_from_center = manifold.center.distance(
						cell->face(f)->vertex(v));
				if (distance_from_center >=  radius + tol) {
					is_inner_face = false;
					break;
				}
			}
			if (is_inner_face) {
				//cout << cell->face(f)->barycenter() << endl;
				cell->face(f)->set_manifold_id(1);
				if (boundary_id > 0){
					cell->face(f)->set_boundary_id(boundary_id);
				}
			}
		}
	}
	mesh.set_manifold(1, manifold);

}

boost::shared_ptr<Mesh<2> > Cylinder2D::makeGrid() {

// this is not the domain by Min and Lee! they had D=1; [-19,50]x[-25,25]
// differences:
// - boundaries are closer (can be fixed with a few additional rectangles)
// - unsymmetric to make a vortex street
// ====> [-19,50]x[-19,25]
//Creation of the principal domain
	const dealii::SphericalManifold<2> manifold1(dealii::Point<2>(0, 0));
	const dealii::SphericalManifold<2> manifold2(dealii::Point<2>(0, 0));
	const dealii::SphericalManifold<2> manifold3(dealii::Point<2>(0, 0));

	dealii::Triangulation<2> mesh;
	dealii::Triangulation<2> rect;
	dealii::Triangulation<2> inner;
	dealii::Triangulation<2> trans;
	dealii::Triangulation<2> tmp;
	dealii::Triangulation<2> merge;
	dealii::GridGenerator::hyper_cube_with_cylindrical_hole(tmp, 5, 19);

// make circular
	make_inner_manifold(tmp, manifold1, 5);
	tmp.refine_global(2);
	tmp.set_manifold(1,manifold1);
// delete coarse cells (Mesh.merge() requires "coarse" triangulations as input)
	dealii::GridGenerator::flatten_triangulation(tmp, trans);
// inner part
	tmp.clear();
	dealii::GridGenerator::hyper_shell(tmp, dealii::Point<2>(0,0), 1, 5, 8);
	tmp.set_all_manifold_ids(1);
	tmp.set_manifold(1, manifold2);
	tmp.refine_global(2);
	tmp.set_manifold(1,manifold2);
	dealii::GridGenerator::flatten_triangulation(tmp, inner);
	dealii::GridGenerator::merge_triangulations(trans, inner, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);
// tail part
	std::vector<unsigned int> repetitions;
	repetitions.push_back(3);
	repetitions.push_back(8);
	dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
			dealii::Point<2>(19, -19), dealii::Point<2>(50, 19));
	dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);
// add unsymmetric part
// rect left
	repetitions.at(0) = 8;
	repetitions.at(1) = 1;
	rect.clear();
	dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
			dealii::Point<2>(-19, 19), dealii::Point<2>(19, 25));
	dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);
// rect right 1
	repetitions.at(0) = 3;
	repetitions.at(1) = 1;
	rect.clear();
	dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
			dealii::Point<2>(19, 19), dealii::Point<2>(50, 25));
	dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);

// set boundary ids
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 0, -19, 1); // left
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 0, 50, 2);
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 1, -19, 3);
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 1, 25, 4);
	make_inner_manifold(mesh, manifold3, 0.5, 0);
	mesh.set_manifold(1,manifold3);

	std::stringstream s;
	s << "/tmp/grid_cylinder.vtk";
	cout << s.str() << endl;
	std::ofstream out(s.str());
	dealii::GridOut grid_out;
	grid_out.write_vtk(mesh, out);
	out.close();

	/*
	 dealii::GridIn<2> gridin;
	 gridin.attach_triangulation(*tria);
	 // !!! When including unv meshes from Salome, the first two blocks have to be deleted manually (compare Original_Mesh_1.unv to Mesh_1.unv)
	 // The mesh file then regularly begins with the lines "-1" and "2411"
	 std::stringstream s;
	 s << getenv("NATRIUM_DIR") << "/src/examples/step-9/salome/Mesh_1.unv";
	 std::ifstream f(s.str().c_str());
	 gridin.read_unv(f);
	 CFDSolverUtilities::mesh_info(*tria, "cylinder.eps");
	 */

	boost::shared_ptr<Mesh<2> > distributed_mesh = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
	distributed_mesh->copy_triangulation(mesh);

// release manifolds
	//tmp.set_manifold(0);

	return distributed_mesh;
}

boost::shared_ptr<BoundaryCollection<2> > Cylinder2D::makeBoundaries(
		double inletVelocity) {

// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = inletVelocity;

	dealii::Tensor<1,2> inletTensor();

	/*boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(4,
					constantVelocity)); */
	boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(1, constantVelocity));
	boundaries->addBoundary(
			boost::make_shared<DoNothingBoundary<2> >(2));
    boundaries->addBoundary(
            boost::make_shared<PeriodicBoundary<2> >(3, 4, 1, getMesh()));

/*	boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(3,
					constantVelocity)); */
	boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(0, zeroVelocity));

// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}


/**
 * @short set initial velocities
 *
 */
double Cylinder2D::InitialU::value(const dealii::Point<2>& x,
			const unsigned int component) const {
		assert(component < 2);
	if (0 == component){

		return m_flow->getCharacteristicVelocity();// * ( 1 + 0.5 * sin(x(1))*cos(5*x(1)+x(0)));

	}
    if (1 == component){
        return m_flow->getCharacteristicVelocity()*0.1 * sin(x(1))*cos(5*x(1)+x(0));// * ( 1 + 0.5 * sin(x(1))*cos(5*x(1)+x(0)));
    }
	return 0;
}

} /* namespace natrium */
