/**
 * @file LidDrivenCavity2D.h
 * @short Lid-driven cavity with three static walls and one moving wall
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef LIDDRIVENCAVIT2D_H_
#define LIDDRIVENCAVIT2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {


/**
 * @short Description of a lid-driven cavity flow.
 */
class LidDrivenCavity2D: public ProblemDescription<2> {
private:
	double topPlateVelocity;
public:

	/// constructor
	LidDrivenCavity2D(double velocity, double viscosity, size_t refinementLevel);

	/// destructor
	virtual ~LidDrivenCavity2D();

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled;
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 * @note Usually, the densities are set to unity.
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

	virtual double getCharacteristicVelocity() const {
		return topPlateVelocity;
	}

private:

	/**
	 * @short create triangulation for lid-driven cavity flow.
	 * @param refinementLevel (denoted as N) The grid will have 2^n*2^n even-sized square cells
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > makeGrid(size_t refinementLevel);

	/**
	 * @short create boundaries for lid-driven cavity flow
	 * @return shared pointer to BoundaryCollection. BoundaryCollection is a container class for
	 *         the specified boundary conditions.
	 * @note All boundary types are inherited of the class Boundary
	 *       Here, we have four wall boundaries. One of them (the upper one)
	 *       has a tangential speed.
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries();

};

} /* namespace natrium */
#endif /* LIDDRIVENCAVIT2D_H_ */
