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

#include "natrium/problemdescription/PeriodicBoundary.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class SteadyPeriodicTestFlow2D: public ProblemDescription<2> {
private:
	size_t m_refinementLevel;
public:
	// class that represents the initial velocity
	class InitialVelocity: public dealii::Function<2> {
	public:
		InitialVelocity(){
		}
		virtual double value(const dealii::Point<2>& ,
				const unsigned int component = 0) const {
			assert(component < 2);
			return 0.1;

		}
	};

	/// constructor
	SteadyPeriodicTestFlow2D(double viscosity, size_t refinementLevel) :
			ProblemDescription<2>(makeGrid(), viscosity, 1), m_refinementLevel(refinementLevel) {
		setInitialU(boost::make_shared<InitialVelocity>());
		/// apply boundary values
		setBoundaries(makeBoundaries());
		// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell


	}

	/// destructor
	virtual ~SteadyPeriodicTestFlow2D() {
	}
	virtual void refine(){
		getMesh()->refine_global(m_refinementLevel);

	}
	virtual void transform(Mesh<2>& mesh){

	}


private:

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid() {
		//Creation of the principal domain
		boost::shared_ptr<Mesh<2> > unitSquare = boost::make_shared<Mesh<2> >(
#ifdef WITH_TRILINOS_MPI
				MPI_COMM_WORLD
#endif
				);
		dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

		// Assign boundary indicators to the faces of the "parent cell"
		Mesh<2>::active_cell_iterator cell = unitSquare->begin_active();
		cell->face(0)->set_all_boundary_ids(0);  // left
		cell->face(1)->set_all_boundary_ids(1);  // right
		cell->face(2)->set_all_boundary_ids(2);  // top
		cell->face(3)->set_all_boundary_ids(3);  // bottom

		return unitSquare;
	}

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries() {

		// make boundary description
		boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
				BoundaryCollection<2> >();
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));

		return boundaries;
	}

};

class UnsteadyPeriodicTestFlow2D: public SteadyPeriodicTestFlow2D {
public:
	// class that represents the initial velocity
	class UnsteadyInitialVelocity: public dealii::Function<2> {
	public:
		UnsteadyInitialVelocity(){
		}
		virtual double value(const dealii::Point<2>& x,
				const unsigned int component = 0) const {
			assert(component < 2);
			if (component == 0) {
				if ((x(1) >= 0.25) and (x(1) < 0.75)) {
					return 0.1;
				} else {
					return -0.1;
				}
			} else {
				return 0.0;
			}
		}
	};

	/// constructor
	UnsteadyPeriodicTestFlow2D(double viscosity, size_t refinementLevel) :
			SteadyPeriodicTestFlow2D(viscosity, refinementLevel) {
		setInitialU(boost::make_shared<UnsteadyInitialVelocity>());
	}
	/// destructor
	virtual ~UnsteadyPeriodicTestFlow2D() {
	}

};

} /* namespace natrium */
#endif /* PeriodicTestFlow2D_H_ */
