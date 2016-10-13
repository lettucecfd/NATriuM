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
	Cylinder2D(double viscosity, double inletVelocity);

	/// destructor
	virtual ~Cylinder2D();

	virtual double getCharacteristicVelocity() const {
		return m_inletVelocity;
	}

	virtual void refine(Mesh<2>& mesh) {
		// remember to make_inner_manifold when implementing global refinement
	}

	virtual void transform(Mesh<2>& mesh) {

	}

private:

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

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialDensities(distributed_vector& initialDensities,
			const vector<dealii::Point<2> >& supportPoints) const;

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<2> >& supportPoints) const;

};

} /* namespace natrium */
#endif /* CYLINDER2D_H_ */
