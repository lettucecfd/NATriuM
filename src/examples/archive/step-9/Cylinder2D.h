/**
 * @file Cylinder2D.h
 * @short Description of the circular cylinder benchmark (Karman vortex street)
 * @date 09.10.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CYLINDER2D_H_
#define CYLINDER2D_H_

#include "deal.II/grid/tria.h"

#include <deal.II/grid/manifold_lib.h>

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/**
 * @short Description of the flow around a circular cylinder (regular channel flow in square domain).
 */
class Cylinder2D: public ProblemDescription<2> {
public:

	/// constructor
	Cylinder2D(double viscosity, double inletVelocity,size_t refinementLevel);

	/// destructor
	virtual ~Cylinder2D();

	virtual double getCharacteristicVelocity() const {
		return m_inletVelocity;
	}

	virtual void refine(Mesh<2>& mesh) {
		// remember to make_inner_manifold when implementing global refinement
		const dealii::SphericalManifold<2> manifold(dealii::Point<2>(0, 0));
		for(size_t i = 0; i<m_refinementLevel;i++) {
            make_inner_manifold(mesh, manifold);
            mesh.refine_global(1);
            // free the manifold
            mesh.set_manifold(1, manifold);
        }
	}

	virtual void transform(Mesh<2>& ) {

	}

	virtual bool isCartesian(){
		return false;
	}

	/**
	 * @short class to describe the the initial solution
	 */
	class InitialU: public dealii::Function<2> {
	private:
		Cylinder2D* m_flow;
	public:
		InitialU(Cylinder2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

    class InitialT: public dealii::Function<2> {
    private:
        Cylinder2D* m_flow;
    public:
        InitialT(Cylinder2D* flow) :
                m_flow(flow) {
        }
        virtual double value(const dealii::Point<2>& x,
                             const unsigned int component = 0) const;
    };

    class InitialRho: public dealii::Function<2> {
    private:
        Cylinder2D* m_flow;
    public:
        InitialRho(Cylinder2D* flow) :
                m_flow(flow) {
        }
        virtual double value(const dealii::Point<2>& x,
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
	static void make_inner_manifold(dealii::Triangulation<2>& mesh,
			const dealii::SphericalManifold<2>& manifold, double radius = 0.5, int boundary_id = -1,
			double tol = 1e-10);

	/**
	 * @short create triangulation for cylinder flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid();

	/**
	 * @short create boundaries for cylinder flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries(
			double inletVelocity);



};

} /* namespace natrium */
#endif /* CYLINDER2D_H_ */
