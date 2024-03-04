/**
 * @file Sphere.h
 * @short Flow around a sphere
 * @date 15.02.2022
 * @author Dominik Wilde, UC San Diego
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "deal.II/grid/tria.h"

#include <deal.II/grid/manifold_lib.h>

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/**
 * @short Description of the flow around a circular cylinder (regular channel flow in square domain).
 */
    class Sphere: public ProblemDescription<3> {
    public:

        /// constructor
        Sphere(double viscosity, double inletVelocity,size_t refinementLevel);

        /// destructor
        virtual ~Sphere();

        virtual double getCharacteristicVelocity() const {
            return m_inletVelocity;
        }

        virtual void refine(Mesh<3>& mesh) {
            // remember to make_inner_manifold when implementing global refinement
            const dealii::SphericalManifold<3> manifold(dealii::Point<3>(0, 0, 0));
                mesh.reset_all_manifolds();
                mesh.set_all_manifold_ids(0);
                make_inner_manifold(mesh, manifold, 0.5, 1);
                dealii::TransfiniteInterpolationManifold<3> transfinite_manifold;
                transfinite_manifold.initialize(mesh);

                mesh.set_manifold(0, transfinite_manifold);
                mesh.set_manifold(1, manifold);

                mesh.refine_global(m_refinementLevel);


        }

        virtual void transform(Mesh<3>& ) {

        }

        virtual bool isCartesian(){
            return false;
        }

        /**
         * @short class to describe the the initial solution
         */
        class InitialU: public dealii::Function<3> {
        private:
            Sphere* m_flow;
        public:
            InitialU(Sphere* flow) :
                    m_flow(flow) {
            }
            virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;
        };

        class InitialT: public dealii::Function<3> {
        private:
            Sphere* m_flow;
        public:
            InitialT(Sphere* flow) :
                    m_flow(flow) {
            }
            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };

        class InitialRho: public dealii::Function<3> {
        private:
            Sphere* m_flow;
        public:
            InitialRho(Sphere* flow) :
                    m_flow(flow) {
            }
            virtual double value(const dealii::Point<3>& x,
                                 const unsigned int component = 0) const;
        };

    private:
        const size_t m_refinementLevel;

        const double m_inletVelocity;

        /**
         * @ short Sets the manifold id and boundary id of the inner manifold, respectively
         * @ short The manifold will contain all faces that have a distance <= radius + tol to the center
         * @param[in] mesh the mesh
         * @param[in] manifold the spherical manifold (that also contains the center)
         * @param[in] radius distance from center to sphere (default: 0.5)
         * @param[in] boundary_id: if >= 0, then all faces in the manifold get this boundary id
         * @param[in] tol tolerance for checking the distance (default: 1e-10)
         */
        static void make_inner_manifold(dealii::Triangulation<3>& mesh,
                                        const dealii::SphericalManifold<3>& manifold, double radius = 0.5, int boundary_id = -1,
                                        double tol = 1e-10);

        /**
         * @short create triangulation for cylinder flow
         * @return shared pointer to a triangulation instance
         */
        boost::shared_ptr<Mesh<3> > makeGrid();

        /**
         * @short create boundaries for cylinder flow
         * @return shared pointer to a vector of boundaries
         * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
         */
        boost::shared_ptr<BoundaryCollection<3> > makeBoundaries(
                double inletVelocity);



    };

} /* namespace natrium */
#endif /* SPHERE_H_ */
