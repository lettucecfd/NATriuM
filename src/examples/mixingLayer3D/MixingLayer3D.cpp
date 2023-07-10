/*
 * MixingLayer3D.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */

#include "MixingLayer3D.h"

#include "deal.II/grid/grid_generator.h"
//#include "deal.II/grid/tria_accessor.h"
//#include "deal.II/grid/tria_iterator.h"
//#include "deal.II/grid/grid_out.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"

float shearlayerthickness = 0.093;

namespace natrium {

    MixingLayer3D::MixingLayer3D(double viscosity,
                                             size_t refinementLevel, double cs) :
            ProblemDescription<3>(makeGrid(), viscosity, 1), m_cs(cs), m_refinementLevel(refinementLevel) {


        /// apply boundary values
        setBoundaries(makeBoundaries());
        // apply analytic solution
        this->setInitialU(boost::make_shared<InitialVelocity>(this));
        this->setInitialRho(boost::make_shared<InitialDensity>(this));
        this->setInitialT(boost::make_shared<InitialTemperature>(this));


    }

    MixingLayer3D::~MixingLayer3D() = default;

    double MixingLayer3D::InitialVelocity::value(const dealii::Point<3>& x, const unsigned int component) const {

        assert(component < 3);
        if (component == 0) {
            return tanh(-x(1)/(2*shearlayerthickness));
        } else {
            return 0.0;
        }
    }

    double MixingLayer3D::InitialDensity::value(const dealii::Point<3>& x, const unsigned int component) const {
        assert(component == 0);
        return 1.0;// + p / (m_flow->m_cs * m_flow->m_cs);
    }

    double MixingLayer3D::InitialTemperature::value(const dealii::Point<3>& x, const unsigned int component) const {
        assert(component == 0);
        return 1.0;
    }

/**
 * @short create triangulation for Compressible Mixing Layer flow
 * @return shared pointer to a triangulation instance
 */
    boost::shared_ptr<Mesh<3> > MixingLayer3D::makeGrid() {
    //Creation of the principal domain

    boost::shared_ptr<Mesh<3> > cube = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
    double lx = 1720 * shearlayerthickness / 2;
    double ly = 387 * shearlayerthickness / 2;
    double lz = 172 * shearlayerthickness / 2;
    dealii::Point<3> corner1(-lx, -ly, -lz);
    dealii::Point<3> corner2(lx, ly, lz);
    std::vector<unsigned int> rep;
    rep.push_back(1);
    rep.push_back(1);
    rep.push_back(1);
    dealii::GridGenerator::subdivided_hyper_rectangle(*cube, rep, corner1, corner2, true);

    // Assign boundary indicators to the faces of the "parent cell"
    /*Mesh<3>::active_cell_iterator cell = cube->begin_active();
    cell->face(0)->set_all_boundary_ids(0);  //
    cell->face(1)->set_all_boundary_ids(1);  //
    cell->face(2)->set_all_boundary_ids(2);  //
    cell->face(3)->set_all_boundary_ids(3);  //
    cell->face(4)->set_all_boundary_ids(4);  //
    cell->face(5)->set_all_boundary_ids(5);  //
*/
    return cube;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3> > MixingLayer3D::makeBoundaries() {

    // make boundary description
    boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
            BoundaryCollection<3> >();

    // velocity vector moving forward
    dealii::Vector<double> plusVector(3);
    plusVector[0]=1.0;
    plusVector[1]=0.0;
    plusVector[2]=0.0;

    // velocity vector moving backward
    dealii::Vector<double> minusVector(3);
    minusVector[0]=-1.0;
    minusVector[1]=0.0;
    minusVector[2]=0.0;

    // set boundaries on top and bottom to move forward / backward
    boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3> >(2, plusVector));
    boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3> >(3, minusVector));

    // set a boundary between 0 and 1, and 4 and 5, with direction 0 (x) and 2 (z), respectively
    boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
    boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(4, 5, 2, getMesh()));

    // Get the triangulation object (which belongs to the parent class).
    boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();
    return boundaries;
}
} /* namespace natrium */

