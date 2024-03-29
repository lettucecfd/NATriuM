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

Cylinder2D::Cylinder2D(double viscosity, double inletVelocity, size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_inletVelocity(
				inletVelocity), m_refinementLevel(refinementLevel) {
	setCharacteristicLength(1.0);

	this->setInitialU(boost::make_shared<InitialU>(this));
    this->setInitialT(boost::make_shared<InitialT>(this));
    this->setInitialRho(boost::make_shared<InitialRho>(this));
    pout << "Initial values set \n";

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
				if (boundary_id >= 0){
					cell->face(f)->set_boundary_id(boundary_id);
					pout << "Face set to boundary id" << boundary_id << ". \n" ;
				}
			}
		}
	}
	mesh.set_manifold(1, manifold);

}

boost::shared_ptr<Mesh<2> > Cylinder2D::makeGrid() {

    const double upper_extent = 6;
    const double total_length = 30;
    const double total_height = 15;

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
	dealii::GridGenerator::hyper_cube_with_cylindrical_hole(tmp, 4, upper_extent);

// make circular
	//make_inner_manifold(tmp, manifold1, 3);
	tmp.refine_global(2);
	tmp.set_manifold(1,manifold1);
// delete coarse cells (Mesh.merge() requires "coarse" triangulations as input)
	dealii::GridGenerator::flatten_triangulation(tmp, trans);
// inner part
	tmp.clear();
	dealii::GridGenerator::hyper_shell(tmp, dealii::Point<2>(0,0), 1, 4, 4);
	tmp.set_all_manifold_ids(1);
	tmp.set_manifold(1, manifold2);
	tmp.refine_global(3);
	tmp.set_manifold(1,manifold2);
	dealii::GridGenerator::flatten_triangulation(tmp, inner);
	dealii::GridGenerator::merge_triangulations(trans, inner, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);



// tail part
	std::vector<unsigned int> repetitions;
	repetitions.push_back(6);
	repetitions.push_back(8);
	dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
			dealii::Point<2>(upper_extent, -upper_extent), dealii::Point<2>(total_length, upper_extent));
	dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);

// add unsymmetric part
//  UPPER rect left
	 repetitions.at(0) = 8;
	repetitions.at(1) = 4;
	rect.clear();
	dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
			dealii::Point<2>(-upper_extent, upper_extent), dealii::Point<2>(upper_extent, total_height));
	dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);

	// UPPER rect right 1
	repetitions.at(0) = 6;
	repetitions.at(1) = 4;
	rect.clear();
	dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
			dealii::Point<2>(upper_extent, upper_extent), dealii::Point<2>(total_length, total_height));
	dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
	mesh.clear();
	mesh.copy_triangulation(merge);


    // LOWER rect left
    repetitions.at(0) = 8;
    repetitions.at(1) = 4;
    rect.clear();
    dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
                                                      dealii::Point<2>(-upper_extent, -upper_extent), dealii::Point<2>(upper_extent, -total_height));
    dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
    mesh.clear();
    mesh.copy_triangulation(merge);

    // LOWER rect right 1
    repetitions.at(0) = 6;
    repetitions.at(1) = 4;
    rect.clear();
    dealii::GridGenerator::subdivided_hyper_rectangle(rect, repetitions,
                                                      dealii::Point<2>(upper_extent, -upper_extent), dealii::Point<2>(total_length, -total_height));
    dealii::GridGenerator::merge_triangulations(mesh, rect, merge);
    mesh.clear();
    mesh.copy_triangulation(merge);



// set boundary ids
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 0, -upper_extent, 1); // left
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 0, total_length, 2);
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 1, -total_height, 3);
	DealIIExtensions::set_boundary_ids_at_hyperplane<2>(mesh, 1, total_height, 4);
	make_inner_manifold(mesh, manifold3, 1.0, 0);
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

	// Inlet
	boundaries->addBoundary(
			boost::make_shared<SLEquilibriumBoundary<2> >(1, constantVelocity, 1.0));
	// Outlet
	boundaries->addBoundary(
			boost::make_shared<DoNothingBoundary<2> >(2));
	// Top & Bottom
    boundaries->addBoundary(
            boost::make_shared<DoNothingBoundary<2> >(3));
    boundaries->addBoundary(
            boost::make_shared<DoNothingBoundary<2> >(4));
    // Obstacle
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
    // use potential flow for initial values
	const double U = m_flow->getCharacteristicVelocity();
	const double R = 0.5;
	const double r = sqrt(x(0)*x(0)+x(1)*x(1));
	assert(r>=0.97);
	const double sin_phi = x(1)/r;
	const double cos_phi = x(0)/r;

	const double u_r =      U * (1 - R*R/(r*r)) * cos_phi;
	const double u_phi = -  U * (1 + R*R/(r*r)) * sin_phi;


    if (0 == component) {
        return u_r * cos_phi -  u_phi * sin_phi;
        }

    if (1 == component) {
        return u_r * sin_phi +  u_phi * cos_phi;
        }
	return 0;
}



double Cylinder2D::InitialT::value(const dealii::Point<2>& x, const unsigned int component) const {

        return 1.0;

}

    double Cylinder2D::InitialRho::value(const dealii::Point<2>& x, const unsigned int component) const {
        // use potential flow for initial values

        const double U = m_flow->getCharacteristicVelocity()*0.5;
        const double R = 0.5;
        const double r = sqrt(x(0)*x(0)+x(1)*x(1));
        assert(r>=0.97);
        const double sin_phi = x(1)/r;
        const double cos_phi = x(0)/r;
        const double cos_2phi = cos_phi*cos_phi - sin_phi*sin_phi;

        const double u_r =      U * (1 - R*R/r*r) * cos_phi;
        const double u_phi = -  U * (1 + R*R/r*r) * sin_phi;

        return 1.0 + 0.5*U*U*(2*R*R/(r*r)*cos_2phi - R*R*R*R/(r*r*r*r));

    }


} /* namespace natrium */

