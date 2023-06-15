/*
 * MixingLayer3D.h
 *
 *  Created on: Sep 18, 2014
 *      Author: dominik
 */

/**
 * @file MixingLayer3D.h
 * @short Description of a simple Periodic Flow (in cubic domain).
 */

#ifndef MixingLayer3D_H_
#define MixingLayer3D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
    class MixingLayer3D: public ProblemDescription<3> {
    public:
        double getCharacteristicVelocity() const {
            return m_uCl;
        }
        /**
         * @short class to describe the x-component of the initial velocity
         * @note other are default (v0=w0=0, rho0=1)
         */
        class InitialVelocity: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            InitialVelocity(MixingLayer3D* flow) : m_flow(flow) { }
            virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };
        class InitialDensity: public dealii::Function<3> {
        private: MixingLayer3D* m_flow;
        public:
            InitialDensity(MixingLayer3D* flow) : m_flow(flow) { }
            virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };

        class InitialTemperature: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            InitialTemperature(MixingLayer3D* flow) :
                    m_flow(flow) { }
            virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };

        class IncompressibleU: public dealii::Function<3> {
        private:
            MixingLayer3D *m_flow;
            InitialVelocity m_initialU;
            double m_increment;
        public:
            //double *m_maxIncUtrp;
            IncompressibleU(MixingLayer3D *flow) :
                    m_flow(flow), m_initialU(flow), m_increment(1e-6) {
                //m_maxIncUtrp = new double (0.0);
            }
//		~IncompressibleU(){
//		        delete m_maxIncUtrp;
//		}

            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };

        class MeanVelocityProfile: public dealii::Function<3> {
        private:
            MixingLayer3D *m_flow;
            IncompressibleU m_initialIncompressibleU;
        public:
            MeanVelocityProfile(MixingLayer3D *flow) :
                    m_flow(flow), m_initialIncompressibleU(flow) {
            }
            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };

        /// constructor
        MixingLayer3D(double viscosity, size_t refinementLevel, double cs = 0.57735026919,
                               std::vector<unsigned int> repetitions, double ReTau = 180.0,
                               double u_cl = 10, double height = 1.0, double length = 6.0,
                               double width = 3.0,
                               bool is_periodic = true, double gridDensity=0.8);

        /// destructor
        virtual ~MixingLayer3D();

        /////////////////////////////////
        // GETTER     // SETTER        //
        /////////////////////////////////

        double getRefinementLevel() const {
            return m_refinementLevel;
        }

        std::vector<unsigned int> getRepetitions() const {
            return m_repetitions;
        }

        double getFrictionReNumber() const {
            return m_ReTau;
        }

        double getCenterLineVelocity() const {
            return m_uCl;
        }

        double getHeight() const {
            return m_height;
        }

        double getLength() const {
            return m_length;
        }

        double getWidth() const {
            return m_width;
        }

        double getMaxUtrp() const {
            return m_maxUtrp;
        }
        double getMaxIncUtrp() const {
            return m_maxIncUtrp;
        }

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
        double m_uCl;
        size_t m_refinementLevel;
        std::vector<unsigned int> m_repetitions;
        double m_ReTau;
        double m_height;
        double m_length;
        double m_width;
        double m_maxUtrp;
        double m_maxIncUtrp;
        double m_gridDensity;

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
