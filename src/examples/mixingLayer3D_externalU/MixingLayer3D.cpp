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
        double k = 1; // waveVectorMagnitude
        double kZero = 23.66 * shearlayerthickness;
        double du = 2;
        double sine = 0.0;

        double freq[3][20] = {{66.2845548 , 57.4228396 ,  0.27347421, 14.29903971, 36.22500418,
                               32.67204667, 90.67090255, 47.18529938, 14.06418876, 70.97004456,
                               62.29178351, 41.45700503, 19.10988557,  8.97074052, 35.46145436,
                               42.82895781, 83.62131221, 82.69185805, 43.59267645, 23.30396557},
                              {48.63136969, 53.30015739, 65.35687411, 32.56024948, 48.86982954,
                               17.95241707, 10.43078963, 63.22671926, 22.08813559, 48.56074986,
                               60.12496519, 16.33120125, 26.78316967, 88.75881632, 71.64466533,
                               52.34201182, 57.25319528, 18.77557742, 89.64339506, 64.39756756},
                              {37.66260713, 86.91397046, 83.82441442, 99.28408371,  5.6199117 ,
                               53.55657792, 51.90739747, 37.04034452, 17.44584865, 72.59120601,
                               2.11968253, 83.72928205, 11.90074664,  1.31334284, 82.25965319,
                               92.9546491 , 89.77514094, 72.85466858, 62.75285836, 92.14498088}};
        double amp[3][20] = {{-0.04874089, -0.40230095,  0.44781972,  0.26948188, -0.30267951,
                               0.15619175, -0.0472235 , -0.15130775,  0.39634864,  0.21417819,
                               0.48157469, -0.02463747,  0.33387625,  0.38073874,  0.16362254,
                               -0.27262168,  0.21187504,  0.42124131, -0.03299644,  0.35060021},
                              {-0.24247737,  0.4571514 ,  0.02742454,  0.16977676,  0.41402119,
                               0.37254541, -0.45054948,  0.03359502, -0.25827299, -0.21240536,
                               0.2884282 ,  0.07797718,  0.47347099,  0.44127   , -0.170844  ,
                               0.40061025, -0.21370938, -0.47579247,  0.28509434, -0.3662098},
                              {-0.23519186, -0.3248561 , -0.46989698,  0.43446305, -0.23997686,
                               0.39800483, -0.19158887, -0.39072899, -0.1267383 , -0.43100461,
                               0.20452294,  0.35672559, -0.19854072,  0.25452   ,  0.4104372 ,
                               0.01309332, -0.24880793, -0.38311154,  0.43897686,  0.20119984}};
        double phas[3][20] = {{-0.25062147,  0.40273954,  0.67196615, -0.12252237,  0.8043513 ,
                               -0.12413927,  0.6193325 ,  0.62263476, -0.64425497,  0.44998552,
                               -0.55313115,  0.02297279, -0.54176666, -0.98065809,  0.36341102,
                               0.9340466 , -0.38262869,  0.1713398 , -0.25984728, -0.18023933},
                              {-0.61019153, -0.31840257,  0.30723735, -0.00192828,  0.57004758,
                               -0.45774561, -0.79907435, -0.79955856,  0.11416249, -0.62429297,
                               0.01899385, -0.3991707 ,  0.26397155, -0.06765534, -0.13509989,
                               -0.11528443,  0.62020873,  0.54140398,  0.18148938, -0.22456433},
                              {-0.23055065,  0.89575861,  0.54250767,  0.80290944, -0.28150918,
                               0.04415374,  0.82023737,  0.05073249,  0.5880679 ,  0.02301412,
                               0.69611798,  0.47873607,  0.22148708,  0.79260772, -0.66682846,
                               -0.00322379, -0.29419377, -0.02807395,  0.90208379,  0.09878267}};

        if (component == 0) {
            for (int j=0; j<=19; j++) {
                sine += (sin(freq[0][j]*x(0) + phas[0][j]) * amp[0][j]) * exp(-pow((x(1)+0.3)/(2*shearlayerthickness),2));
            }
            return du / 2 * tanh(-x(1)/(2*shearlayerthickness)) + sine; // + u_rand[0]; // holger: u1 = ("U/2) tanh(−x2/δθ (0))
        } else {
            for (int j=0; j<=19; j++) {
                sine += (sin(freq[component][j]*x(0) + phas[component][j]) * amp[component][j]) * exp(-pow((x(1)+0.3)/(2*shearlayerthickness),2));
            }
            return sine;
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

