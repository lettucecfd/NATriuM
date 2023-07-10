/*
 * MixingLayer3D.h
 *
 *  Created on: Sep 18, 2014
 *      Author: dominik
 */

#ifndef MixingLayer3D_solenoidal_H_
#define MixingLayer3D_solenoidal_H_

/**
 * @file MixingLayer3D_solenoidal.h
 * @short Description of a simple Periodic Flow (in cubic domain).
 */

#include "deal.II/grid/tria.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
    class MixingLayer3D: public ProblemDescription<3> {
    public:

        /**
         * @short class to describe the x-component of the initial velocity
         * @note other are default (v0=w0=0, rho0=1)
         */
        class InitialVelocity: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            explicit InitialVelocity(MixingLayer3D* flow) : m_flow(flow) { }
            double value(const dealii::Point<3>& x, unsigned int component) const override;
        };
        void setInitialPsi(boost::shared_ptr<dealii::Function<dim> > ini_psi) {
            m_initialPsi = ini_psi;
        }
        class InitialPsi: public dealii::Tensor<1,3> {
        private:
            MixingLayer3D* m_initialPsi;
        };
        class get_rotation_matrix: public dealii::Tensor<1,3> {
        private:
            MixingLayer3D* m_flow;
        public:
            explicit get_rotation_matrix(MixingLayer3D* flow) : m_flow(flow) { }
            static double value(const dealii::Point<3>& x, unsigned int component) ;
        };
//        dealii::Tensor<1, 3> get_rotation_matrix(const vector <dealii::Tensor<1, 3>> &grad_u);
        class InitialDensity: public dealii::Function<3> {
        private: MixingLayer3D* m_flow;
        public:
            explicit InitialDensity(MixingLayer3D* flow) : m_flow(flow) { }
            double value(const dealii::Point<3>& x, unsigned int component) const override;
        };
        class InitialTemperature: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            explicit InitialTemperature(MixingLayer3D* flow) : m_flow(flow) { }
            double value(const dealii::Point<3>& x, unsigned int component) const override;
        };

        /// constructor
        MixingLayer3D(double viscosity, size_t refinementLevel, double cs = 0.57735026919);

        /// destructor
        ~MixingLayer3D() override;

    private:
//            double RandomVelocity(const dealii::Point<3>& x, unsigned int component);

        void refine(Mesh<3>& mesh) override {
            // Refine grid
            mesh.refine_global(m_refinementLevel);
        }
        void transform(Mesh<3>&) override { }
        bool isCartesian() override { return true; }

    private:
        /// speed of sound
        double m_cs;
        size_t m_refinementLevel;
        /**
         * @short create triangulation for couette flow
         * @return shared pointer to a triangulation instance
         */
        boost::shared_ptr<Mesh<3> > makeGrid();
        /**
         * @short create boundaries for couette flow
         * @return shared pointer to a vector of boundaries
         * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
         */
        boost::shared_ptr<BoundaryCollection<3> > makeBoundaries();

//        void randf_2(int idum, int &iy, vector<int> &iv, double &ran1, int &iseed);

//        void setInitialSines(vector<float> sines);
    };
} /* namespace natrium */

#endif /* MixingLayer3D_solenoidal_H_ */
