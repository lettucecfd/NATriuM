/**
 * @file PeriodicTestFlow2D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PERIODICTESTFLOW2D_H_
#define PERIODICTESTFLOW2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "problemdescription/PeriodicBoundary.h"
#include "problemdescription/ProblemDescription.h"
#include "utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class SteadyPeriodicTestFlow2D: public ProblemDescription<2> {
public:
	/// constructor
	SteadyPeriodicTestFlow2D(double viscosity, size_t refinementLevel) :
			ProblemDescription<2>(makeGrid(refinementLevel), viscosity) {

		/// apply boundary values
		setBoundaries(makeBoundaries());
	}

	/// destructor
	virtual ~SteadyPeriodicTestFlow2D() {
	}

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialDensities(distributed_vector& initialDensities,
			vector<dealii::Point<2> >& supportPoints) const {
		for (size_t i = 0; i < initialDensities.size(); i++) {
			initialDensities(i) = 1.0;
		}
	}

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			vector<dealii::Point<2> >& supportPoints) const {
		assert(
				initialVelocities.at(0).size()
						== initialVelocities.at(1).size());
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
			initialVelocities.at(0)(i) = 0.1;
			initialVelocities.at(1)(i) = 0.1;
		}
	}

private:

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > makeGrid(size_t refinementLevel) {
		//Creation of the principal domain
		shared_ptr<Triangulation<2> > unitSquare =
				make_shared<Triangulation<2> >();
		dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

		// Assign boundary indicators to the faces of the "parent cell"
		Triangulation<2>::active_cell_iterator cell =
				unitSquare->begin_active();
		cell->face(0)->set_all_boundary_indicators(0);  // left
		cell->face(1)->set_all_boundary_indicators(1);  // right
		cell->face(2)->set_all_boundary_indicators(2);  // top
		cell->face(3)->set_all_boundary_indicators(3);  // bottom

		// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
		unitSquare->refine_global(refinementLevel);

		return unitSquare;
	}

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries() {

		// make boundary description
		shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
				BoundaryCollection<2> >();
		boundaries->addBoundary(
				make_shared<PeriodicBoundary<2> >(0, 1, getTriangulation()));
		boundaries->addBoundary(
				make_shared<PeriodicBoundary<2> >(2, 3, getTriangulation()));

		// Get the triangulation object (which belongs to the parent class).
		shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();

		return boundaries;
	}

};



class UnsteadyPeriodicTestFlow2D: public SteadyPeriodicTestFlow2D{
public:
	/// constructor
	UnsteadyPeriodicTestFlow2D(double viscosity, size_t refinementLevel) :
		SteadyPeriodicTestFlow2D(viscosity, refinementLevel) {
	}
	/// destructor
	virtual ~UnsteadyPeriodicTestFlow2D() {
	}

	virtual void applyInitialDensities(distributed_vector& initialDensities,
			vector<dealii::Point<2> >& supportPoints) const {
		for (size_t i = 0; i < initialDensities.size(); i++) {
			initialDensities(i) = 1.0;
		}
	}

	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			vector<dealii::Point<2> >& supportPoints) const {
		assert(
				initialVelocities.at(0).size()
						== initialVelocities.at(1).size());
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
			if ((supportPoints.at(i)(1) > 0.25) and (supportPoints.at(i)(1) < 0.75)){
				initialVelocities.at(0)(i) = 0.1;
			} else {
				initialVelocities.at(0)(i) = -0.1;
			}
			initialVelocities.at(1)(i) = 0.0;
		}
	}
};

} /* namespace natrium */
#endif /* PeriodicTestFlow2D_H_ */
