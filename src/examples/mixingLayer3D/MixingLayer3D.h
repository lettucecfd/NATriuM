/*
 * MixingLayer3D.h
 *
 *  Created on: Sep 18, 2014
 *      Author: dominik
 */

#ifndef MixingLayer3D_H_
#define MixingLayer3D_H_

/**
 * @file MixingLayer3D.h
 * @short Description of a simple Periodic Flow (in cubic domain).
 */

#include "deal.II/grid/tria.h"
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"
#include "deal.II/grid/grid_out.h"
#include <math.h>

#include <deal.II/lac/la_parallel_vector.h>

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
                vector<double> xvec;
                vector<double> yvec;
                vector<double> zvec;
//                DealIIExtensions::LinearAlgebra::distributed::Vector<
//                    DealIIExtensions::LinearAlgebra::distributed::Vector<double>> randomPsi;
                vector< vector< vector< vector<double> > > > randomPsi;
                vector< vector< vector< vector<double> > > > curlOfPsi;
                double minx, miny, minz;
                double maxx, maxy, maxz;
//                double dx, dy, dz;
                double lx, ly, lz;
                int nx, ny, nz;
                int kxmax, kymax, kzmax;
                bool m_print, m_recalculate;
                double InterpolateVelocities(double xq, double yq, double zq, const unsigned int dim) const;
            public:
                InitialVelocity(MixingLayer3D *flow, bool print, bool recalculate, string dirName);
                double value(const dealii::Point<3>& x, const unsigned int component = 0) const override;
                vector<vector<vector<std::complex<double>>>> Fourier3D(const vector<vector<vector<double>>> &in) const;
                vector<vector<vector<double>>> InverseFourier3D(const vector<vector<vector<std::complex<double>>>> &in) const;
        };
        class InitialDensity: public dealii::Function<3> {
            private:
                MixingLayer3D* m_flow;
            public:
                explicit InitialDensity(MixingLayer3D* flow) : m_flow(flow) { }
                virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };
        class InitialTemperature: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            InitialTemperature(MixingLayer3D* flow) : m_flow(flow) { }
            virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };
        bool m_squash;

        /// constructor
        MixingLayer3D(double viscosity, size_t refinementLevel, bool squash, bool print, bool recalculate,
                      string dirName, double U = 1.);
        /// destructor
        virtual ~MixingLayer3D();

        virtual void refine(Mesh<3>& mesh) {
            mesh.refine_global(m_refinementLevel);
//            if (m_squash) {
//                dealii::GridTools::transform(UnstructuredGridFunc(), mesh);
//                std::ofstream out("grid-mixinglayer.eps");
//                dealii::GridOut grid_out;
//                grid_out.write_eps(mesh, out);
//            }
        }
        struct UnstructuredGridFunc {
//            double sigma = 0.1;
            double ymax = 387 * 0.093 / 2;
            double k = 1;
            double trans(double y) const {
//                return y * erf((y-ymax)/sqrt(2*sigma*sigma));
//                return yrel*sqrt(1-k/1.05)
                double yrel = y/ymax;
                return pow(yrel, 3) * ymax; // * abs(yrel)
//                return y;
//                return tanh((y/ymax)/(2 * 0.093));
            }
            dealii::Point<3> operator() (const dealii::Point<3> &in) const {
                return dealii::Point<3>(in(0), trans(in(1)), in(2));
            }
        };
        virtual void transform(Mesh<3>&) {}
        virtual bool isCartesian() {return true;}

    private:
        /// speed of sound
        double m_U;
        double lx, ly, lz;
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
//    protected:
    };

} /* namespace natrium */

#endif /* MixingLayer3D_H_ */
