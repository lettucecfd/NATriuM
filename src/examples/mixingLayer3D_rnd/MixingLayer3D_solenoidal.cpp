/*
 * MixingLayer3D.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */

#include "MixingLayer3D_solenoidal.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/grid_out.h"
#include <deal.II/physics/transformations.h>
#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"
#include <random>
#include<list>

float shearlayerthickness = 0.093; // TODO



namespace natrium {

    MixingLayer3D::MixingLayer3D(double viscosity, size_t refinementLevel, double cs) :
            ProblemDescription<3>(makeGrid(), viscosity, 1), m_cs(cs), m_refinementLevel(refinementLevel) {
        /// apply boundary values
        setBoundaries(makeBoundaries());
        this->setInitialPsi(boost::make_shared<InitialPsi>(this));
//        setInitialPsi(boost::make_shared<InitialPsi>(this))
//        this->setInitialSines(boost::make_shared<m_initial_sines>(this));
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
        // initialize velocities
        double x0 = x(0);
        double x1 = x(1);
        double x2 = x(2);

        if (component == 0) {
            return du / 2 * tanh(-x(1)/(2*shearlayerthickness)) + MixingLayer3D::get_rotation_matrix::value(x, component); // + u_rand[0]; // holger: u1 = ("U/2) tanh(−x2/δθ (0))
        } else {
            return MixingLayer3D::get_rotation_matrix::value(x, component);
        }
    }

    boost::make_shared< dealii::Tensor<1,3> > MixingLayer3D::InitialPsi(boost::make_shared< Mesh<3> > domain) {
        return;
    }

    double MixingLayer3D::get_rotation_matrix::value(const dealii::Point<3> &x, unsigned int component) {
        return this->m_curl(x(component));
    }

    boost::make_shared< dealii::Tensor<1,3> > MixingLayer3D::get_rotation_matrix(const std::vector<dealii::Tensor<1, 3>> &grad_u) {
        const dealii::Tensor<1, 3> curl({
                                                grad_u[2][1] - grad_u[1][2],
                                                grad_u[0][2] - grad_u[2][0],
                                                grad_u[1][0] - grad_u[0][1]});
        const double tan_angle = std::sqrt(curl * curl);
        const double angle     = std::atan(tan_angle);
        if (std::abs(angle) < 1e-9)
        {
            static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            static const dealii::Tensor<2, 3> rot(rotation);
            return rot;
        }
        const dealii::Tensor<1, 3> axis = curl / tan_angle;
        return curl; //Physics::Transformations::Rotations::rotation_matrix_3d(axis, -angle);
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

        boost::shared_ptr< Mesh<3> > domain = boost::make_shared< Mesh<3> >(MPI_COMM_WORLD);
        double lx = 1720 * shearlayerthickness / 2;
        double ly = 387 * shearlayerthickness / 2;
        double lz = 172 * shearlayerthickness / 2;
        dealii::Point<3> corner1(-lx, -ly, -lz);
        dealii::Point<3> corner2(lx, ly, lz);
        std::vector<unsigned int> rep;
        rep.push_back(1);
        rep.push_back(1);
        rep.push_back(1);
        dealii::GridGenerator::subdivided_hyper_rectangle(*domain, rep, corner1, corner2, true);

        // Assign boundary indicators to the faces of the "parent cell"
        /*Mesh<3>::active_cell_iterator cell = cube->begin_active();
        cell->face(0)->set_all_boundary_ids(0);  //
        cell->face(1)->set_all_boundary_ids(1);  //
        cell->face(2)->set_all_boundary_ids(2);  //
        cell->face(3)->set_all_boundary_ids(3);  //
        cell->face(4)->set_all_boundary_ids(4);  //
        cell->face(5)->set_all_boundary_ids(5);  //
    */
        return domain;
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
//        boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();
        return boundaries;
    }
//    void MixingLayer3D::InitialVectorPotential()

//    inline void MixingLayer3D::setInitialSines(vector<float> sines){
//        std::random_device rd;  // Will be used to obtain a seed for the random number engine
//        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//        std::uniform_real_distribution<> amp(-0.5, 0.5);
//        std::uniform_real_distribution<> freq(0., 100.0);
//        std::uniform_real_distribution<> phase(-1.0, 1.0);
//
////        std::list<void> sine(3, 0.0);
////        for (int h=0; h<=2; h++){
////            m_initial_sines[h] = 0;
////        }
//        vector< vector<int> > amplitudes;
//        for (int h=0; h<=2; h++){
//            vector<int> amplitudes[h] = 0.0;
//        }
//        for (int j=0; j<=20; j++) {
//            amplitudes[]
//            for (int h=0; h<=2; h++){
//                sine[h] += sin(freq(gen)*sines[h] + phase(gen)) * amp(gen);
//                sine[h] *= exp(-pow((sines[1]+0.3)/(2*shearlayerthickness),2));
//            }
//        }
//        for (int h=0; h<=2; h++){
//            m_initial_sines[h] += sine[h];
//        }
//    }

//    inline void MixingLayer3D::randf_2(int idum, int &iy, vector<int> &iv, double &ran1, int &iseed){
//
//        int			ia		= 16807;
//        double		im 		= 2147483647;
//        double		am		= 1/im;
//        int 		iq		= 127773;
//        int 		ir		= 2836;
//        int 		ntab	= 32;
//        double 		ndiv	= 1+(im-1)/ntab;
//        double 		eps		= 1.2e-7;
//        double		rnmx	= 1-eps;
//
//        int 			j, k;
//
////initial iseed (idum) is negative
//        if (idum <= 0 || iy == 0){
//            idum = std::max(-idum, 1);
//            for (int j = ntab+8; j >=1; --j){
//                k = floor(idum/iq);
//                idum = ia*(idum-k*iq)-ir*k;
//                if (idum < 0){
//                    idum = idum+im;
//                }
//                if (j <= ntab){
//                    iv[j-1] = idum;
//                }
//            }
//            iy = iv[0];
//        }
//
//        k 		= floor(idum/iq);
//        idum 	= ia*(idum-k*iq)-ir*k;
//
//        if  (idum <= 0){
//            idum = idum+im;
//        }
//
//        j 		= floor(iy/ndiv);
//        iy 		= iv[j];
//        iv[j]	= idum;
//        iseed	= idum;
//        ran1	= std::min(am*iy,rnmx);
//    }
//
//    double MixingLayer3D::m_initial_sines::value(const dealii::Point<3> &x, unsigned int component) const {
//        return Function::value(x, component);
//    }
} /* namespace natrium */


//        // setup randomized vector potential field
//        vector<int> psi (3);
//        // psi = {rand() / (RAND_MAX), rand() / (RAND_MAX), rand() / (RAND_MAX)};
//        for (int & i : psi) {
//            psi[i] = ((double) rand() / (RAND_MAX)) * exp(-2*k/kZero) * exp(-pow((x(1)+0.3)/(2*shearlayerthickness),2));
//        }
//         x = x(0), y = x(1), z = x(2)
//         psi_x = psi[0] // R^3
//        for (int dim=0; 3; dim++) {
//
//        }
//        // dpsix_dy = { psi[0][:,i,:]-psi[0][:,i+1,:] for i in psi[0].shape(1) }
//        std::vector<std::vector<double>>> dpsi (3);
//        dpsi[0] = {dpsix_dx, dpsix_dy, dpsix_dz};
//        dpsi[1] = {dpsiy_dx, dpsiy_dy, dpsiy_dz};
//        dpsi[2] = {dpsiz_dx, dpsiz_dy, dpsiz_dz};
//        // calculate rotation of vector potentail field
//        vector<int> u_rand (3);
//        for (int & i : u_rand) {
//            u_rand[i] = dpsi[i-1,i+1] - dpsi[i+1,i-1];
//        }
//        double centering = exp(-2*k/kZero) * exp(-pow((x(1))/(2*shearlayerthickness),2));
//        vector<int> u_rand (3);
//        for (int & i : u_rand) {
//            u_rand[i] = ((double) rand() / (RAND_MAX)) * centering;
//        }
//        std::random_device rd;  // Will be used to obtain a seed for the random number engine
//        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//        std::uniform_real_distribution<> amp(-1.0, 1.0);
//        std::uniform_real_distribution<> freq(0., 100.0);
//        std::uniform_real_distribution<> phase(-1.0, 1.0);
//
//        vector<int> u_rand (3);
//        for (int & i : u_rand) {
//            u_rand[i] = 0.;
//            for (int j=0; j<=10; j++) {
//                double sine;
//                for (int h=0; h<=2; h++){
//                    sine += sin(freq(gen)*x(h) + phase(gen));
//                }
//                sine *= exp(-pow((x(1)+0.3)/(2*shearlayerthickness),2)) * amp(gen);
//                u_rand[i] += sine;
//            }
//        }