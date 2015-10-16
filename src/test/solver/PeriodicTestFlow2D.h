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
public:
	// class that represents the initial velocity
	class InitialVelocity: public dealii::Function<2> {
	private:
		SteadyPeriodicTestFlow2D* m_flow;
	public:
		InitialVelocity(SteadyPeriodicTestFlow2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x,
				const unsigned int component = 0) const {
			assert(component < 2);
			return 0.1;

		}
	};

	/// constructor
	SteadyPeriodicTestFlow2D(double viscosity, size_t refinementLevel) :
			ProblemDescription<2>(makeGrid(), viscosity, 1) {
		setInitialU(make_shared<InitialVelocity>(this));
		/// apply boundary values
		setBoundaries(makeBoundaries());
		// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
		getMesh()->refine_global(refinementLevel);

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
			const map<dealii::types::global_dof_index, dealii::Point<2> >& ) const {
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
			const map<dealii::types::global_dof_index, dealii::Point<2> >& ) const {
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
	shared_ptr<Mesh<2> > makeGrid() {
		//Creation of the principal domain
		shared_ptr<Mesh<2> > unitSquare = make_shared<Mesh<2> >(
#ifdef WITH_TRILINOS_MPI
				MPI_COMM_WORLD
#endif
				);
		dealii::GridGenerator::hyper_cube(*unitSquare, 0, 1);

		// Assign boundary indicators to the faces of the "parent cell"
		Mesh<2>::active_cell_iterator cell = unitSquare->begin_active();
		cell->face(0)->set_all_boundary_indicators(0);  // left
		cell->face(1)->set_all_boundary_indicators(1);  // right
		cell->face(2)->set_all_boundary_indicators(2);  // top
		cell->face(3)->set_all_boundary_indicators(3);  // bottom

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
				make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
		boundaries->addBoundary(
				make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));

		// Get the triangulation object (which belongs to the parent class).
		shared_ptr<Mesh<2> > tria_pointer = getMesh();

		return boundaries;
	}

};

class UnsteadyPeriodicTestFlow2D: public SteadyPeriodicTestFlow2D {
public:
	// class that represents the initial velocity
	class UnsteadyInitialVelocity: public dealii::Function<2> {
	private:
		UnsteadyPeriodicTestFlow2D* m_flow;
	public:
		UnsteadyInitialVelocity(UnsteadyPeriodicTestFlow2D* flow) :
				m_flow(flow) {
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
		setInitialU(make_shared<UnsteadyInitialVelocity>(this));
	}
	/// destructor
	virtual ~UnsteadyPeriodicTestFlow2D() {
	}

};

} /* namespace natrium */
#endif /* PeriodicTestFlow2D_H_ */
