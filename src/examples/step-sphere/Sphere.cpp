/**
 * @file Sphere.cpp
 * @short Flow around a sphere
 * @date 15.02.2022
 * @author Dominik Wilde, UC San Diego
 */

#include "Sphere.h"

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

    Sphere::Sphere(double viscosity, double inletVelocity, size_t refinementLevel) :
            ProblemDescription<3>(makeGrid(), viscosity, 1.0), m_inletVelocity(
            inletVelocity), m_refinementLevel(refinementLevel) {
        setCharacteristicLength(1.0);

        this->setInitialU(boost::make_shared<InitialU>(this));
        this->setInitialT(boost::make_shared<InitialT>(this));
        this->setInitialRho(boost::make_shared<InitialRho>(this));
        pout << "Initial values set \n";

        /// apply boundary values
        setBoundaries(makeBoundaries(inletVelocity));
    }

    Sphere::~Sphere() {
    }

    void Sphere::make_inner_manifold(dealii::Triangulation<3>& mesh,
                                         const dealii::SphericalManifold<3>& manifold, double radius, int boundary_id,
                                         double tol) {

        dealii::Triangulation<3>::active_cell_iterator cell = mesh.begin_active(),
                endc = mesh.end();
        for (; cell != endc; ++cell) {
            for (unsigned int f = 0; f < dealii::GeometryInfo<3>::faces_per_cell;
                 ++f) {
                bool is_inner_face = true;
                for (unsigned int v = 0;
                     v < dealii::GeometryInfo<3>::vertices_per_face; ++v) {
                    const double distance_from_center = manifold.center.distance(
                            cell->face(f)->vertex(v));
                    if (std::abs(cell->face(f)->vertex(v).norm_square() - 0.25) > 1e-12) {
                        is_inner_face = false;
                        break;
                    }
                }
                if (is_inner_face) {
                    cell->face(f)->set_all_manifold_ids(1);
                    if (boundary_id >= 0){
                        cell->face(f)->set_boundary_id(boundary_id);
                        pout << "Face set to boundary id" << boundary_id << ". \n" ;
                    }
                }
                else {
                    cell->face(f)->set_all_manifold_ids(0);
                }
            }
        }

    }

    boost::shared_ptr<Mesh<3> > Sphere::makeGrid() {
        const double cube_half_width = std::sqrt(3);
        dealii::Triangulation<3> mesh;

        dealii::Triangulation<3> tria_outer;
        dealii::GridGenerator::hyper_shell(
                mesh, dealii::Point<3>(), 0.5, cube_half_width, 6); // 6 because of the cubic shape!
        mesh.reset_all_manifolds();
        mesh.set_all_manifold_ids(0);
        const dealii::SphericalManifold<3> manifold1(dealii::Point<3>(0, 0, 0));
        make_inner_manifold(mesh, manifold1, 0.5,1);


        mesh.set_manifold(1, manifold1);

        dealii::TransfiniteInterpolationManifold<3> transfinite_manifold;
        transfinite_manifold.initialize(mesh);
        mesh.set_manifold(0, transfinite_manifold);




        const int inlet_x = -2;
        const int outlet_x = 8;
        const int width = 5;


        std::vector<unsigned int> repetitions;


        dealii::Triangulation<3> merge;


// make adjacent regions

        std::vector<std::pair<int, int> > lengths = { {inlet_x, -1}, {-1,1}, {1,3},  {3,outlet_x} };
        std::vector<std::pair<int, int> > sizes = { {-width, -1}, {-1,1}, {1,width} };
        for(auto i : lengths){
            for(auto j : sizes){
                for(auto k : sizes) {
                    if(i.first == -1 && j.first == -1 && k.first == -1)
                        continue;
                    dealii::Triangulation<3> new_mesh;
                    dealii::GridGenerator::hyper_rectangle(new_mesh,
                                                           dealii::Point<3>(i.first, j.first, k.first),
                                                           dealii::Point<3>(i.second, j.second, k.second));
                    new_mesh.reset_all_manifolds();
                    new_mesh.set_all_manifold_ids(2);

                    dealii::GridGenerator::merge_triangulations(new_mesh,mesh, mesh, 1.0e-12);



                }
            }
        }








// set boundary ids
        DealIIExtensions::set_boundary_ids_at_hyperplane<3>(mesh, 0, inlet_x, 2); // left
        DealIIExtensions::set_boundary_ids_at_hyperplane<3>(mesh, 0, outlet_x, 3);
        DealIIExtensions::set_boundary_ids_at_hyperplane<3>(mesh, 1, -width, 4);
        DealIIExtensions::set_boundary_ids_at_hyperplane<3>(mesh, 1, width, 5);
        DealIIExtensions::set_boundary_ids_at_hyperplane<3>(mesh, 2, -width, 6);
        DealIIExtensions::set_boundary_ids_at_hyperplane<3>(mesh, 2, width, 7);





        boost::shared_ptr<Mesh<3> > distributed_mesh = boost::make_shared<Mesh<3> >(
                MPI_COMM_WORLD);
        distributed_mesh->copy_triangulation(mesh);
        make_inner_manifold(*distributed_mesh, manifold1, 0.5,1);
        mesh.set_manifold(1, manifold1);


        std::stringstream s;
        s << "./grid_sphere.vtk";
        cout << s.str() << endl;
        std::ofstream out(s.str());
        dealii::GridOut grid_out;
        grid_out.write_vtk(mesh, out);
        out.close();

        return distributed_mesh;
    }

    boost::shared_ptr<BoundaryCollection<3> > Sphere::makeBoundaries(
            double inletVelocity) {

// make boundary description
        boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
                                                               BoundaryCollection<3> >();
        numeric_vector zeroVelocity(3);
        numeric_vector constantVelocity(3);
        constantVelocity(0) = inletVelocity;

        // Inlet
        boundaries->addBoundary(
                boost::make_shared<SLEquilibriumBoundary<3> >(2, constantVelocity, 1.0));
        // Outlet
        boundaries->addBoundary(
                boost::make_shared<SLEquilibriumBoundary<3> >(3, constantVelocity, 1.0));
        // Top & Bottom
        boundaries->addBoundary(
                boost::make_shared<SLEquilibriumBoundary<3> >(4, constantVelocity, 1.0));
        boundaries->addBoundary(
                boost::make_shared<SLEquilibriumBoundary<3> >(5, constantVelocity, 1.0));
        boundaries->addBoundary(
                boost::make_shared<SLEquilibriumBoundary<3> >(6, constantVelocity, 1.0));
        boundaries->addBoundary(
                boost::make_shared<DoNothingBoundary<3> >(7));
        // Obstacle
        boundaries->addBoundary(
                boost::make_shared<VelocityNeqBounceBack<3> >(1, zeroVelocity));

// Get the triangulation object (which belongs to the parent class).
        //boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();

        return boundaries;
    }


/**
 * @short set initial velocities
 *
 */
    double Sphere::InitialU::value(const dealii::Point<3>& x,
                                       const unsigned int component) const {
        assert(component < 4);
        // use potential flow for initial values
        const double U = m_flow->getCharacteristicVelocity();
        const double R = 0.5;
        const double r = sqrt(x(0)*x(0)+x(1)*x(1)+x(2)*x(2));
        const double yz_distance = sqrt(x(1)*x(1)+x(2)*x(2));
        const double sin_phi = yz_distance/r;
        const double cos_phi = x(0)/r;

        const double u_r =      U * (1 - R*R/(r*r)) * cos_phi;
        const double u_phi = -  U * (1 + R*R/(r*r)) * sin_phi;


        if (0 == component) {
            //return u_r * cos_phi -  u_phi * sin_phi;
            return m_flow->getCharacteristicVelocity()/1.7*1.3*(1-sin(x(0)*x(1))*0.1);;
            //U*(1-3*R*R*R*x(0)*x(0)/(2*r*r*r*r*r)+(R*R*R)/(2*r*r*r));

        }


        if (1 == component) {
           // return 0.0;// (u_r * sin_phi +  u_phi * cos_phi) * x(1)/yz_distance;
           return 0.3*sin(x(1)*x(2))*cos(x(0));//-U*3*(R*R*R)/(2*r*r*r*r*r)*x(0)*x(1);
        }

        if (2 == component) {
            //return 0.0; //(u_r * sin_phi +  u_phi * cos_phi) * x(2)/yz_distance;
            return 0.3*sin(x(1)+0.33*x(2)+0.12)*cos(x(0)+0.33);//-U*3*(R*R*R)/(2*r*r*r*r*r)*x(0)*x(2);

        }
        return 0;
    }



    double Sphere::InitialT::value(const dealii::Point<3>& x, const unsigned int component) const {
        const double U = m_flow->getCharacteristicVelocity();
        const double R = 0.5;
        const double r = sqrt(x(0)*x(0)+x(1)*x(1)+x(2)*x(2));
        const double ux = U*(1-3*R*R*R*x(0)*x(0)/(2*r*r*r*r*r)+(R*R*R)/(2*r*r*r));
        const double uy = -U*3*(R*R*R)/(2*r*r*r*r*r)*x(0)*x(1);
        const double uz = -U*3*(R*R*R)/(2*r*r*r*r*r)*x(0)*x(2);
        const double speed2 = ux*ux+uy*uy+uz*uz;
        return 1.0;//1.5 - speed2/2;

    }

    double Sphere::InitialRho::value(const dealii::Point<3>& x, const unsigned int component) const {
        // use potential flow for initial values
        const double U = m_flow->getCharacteristicVelocity();
        const double R = 0.5;
        const double r = sqrt(x(0)*x(0)+x(1)*x(1)+x(2)*x(2));
        const double ux = U*(1-3*R*R*R*x(0)*x(0)/(2*r*r*r*r*r)+(R*R*R)/(2*r*r*r));
        const double uy = -U*3*(R*R*R)/(2*r*r*r*r*r)*x(0)*x(1);
        const double uz = -U*3*(R*R*R)/(2*r*r*r*r*r)*x(0)*x(2);
        const double speed2 = ux*ux+uy*uy+uz*uz;
        return 1.0;//1.5 - speed2/2;

    }


} /* namespace natrium */

