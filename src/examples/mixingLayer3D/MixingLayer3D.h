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
            InitialVelocity(MixingLayer3D* flow) :
                    m_flow(flow) {
            }
            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };
        class InitialDensity: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            InitialDensity(MixingLayer3D* flow) :
                    m_flow(flow) {
            }
            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };

        class InitialTemperature: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            InitialTemperature(MixingLayer3D* flow) :
                    m_flow(flow) {
            }
            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };

        /// constructor
        MixingLayer3D(double viscosity, size_t refinementLevel, double cs =
        0.57735026919);

        /// destructor
        virtual ~MixingLayer3D();

        virtual void refine(Mesh<3>& mesh) {
            // Refine grid
            mesh.refine_global(m_refinementLevel);
        }
        virtual void transform(Mesh<3>&) {

        }
        virtual bool isCartesian() {
            return true;
        }
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

    };

} /* namespace natrium */

#endif /* MixingLayer3D_H_ */
