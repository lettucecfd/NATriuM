/**
 * @file TaylorGreenTest2D.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TaylorGreenTest2D_H_
#define TaylorGreenTest2D_H_

#include "deal.II/grid/tria.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"


using namespace natrium;


/**
 *  @short Description of a Taylor-Green vortex, a benchmark with only periodic boundaries
 */
class TaylorGreenTest2D: public ProblemDescription<2> {
public:

	/// constructor
	TaylorGreenTest2D(double viscosity, size_t refinementLevel) :
			ProblemDescription<2>(makeGrid(), viscosity, 1), m_refinementLevel(refinementLevel) {

		/// apply boundary values
		setBoundaries(makeBoundaries());

	}

	/// destructor
	virtual ~TaylorGreenTest2D() {
	}
	;

	static void getAnalyticSolution(double time, distributed_vector& analyticSolution1,
			distributed_vector& analyticSolution2,
			const vector<dealii::Point<2> >& supportPoints,
			const TaylorGreenTest2D& tgVortex) {
		assert(analyticSolution1.size() == supportPoints.size());
		assert(analyticSolution2.size() == supportPoints.size());
		assert(supportPoints.size() > 0);

		for (size_t i = 0; i < supportPoints.size(); i++) {
			analyticSolution1(i) = tgVortex.analyticVelocity1(supportPoints.at(i),
					time);
			analyticSolution2(i) = tgVortex.analyticVelocity2(supportPoints.at(i),
					time);
		}
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
			const map<dealii::types::global_dof_index, dealii::Point<2> >& supportPoints) const {
		assert(
				initialVelocities.at(0).size()
						== initialVelocities.at(1).size());
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
			initialVelocities.at(0)(i) = analyticVelocity1(supportPoints.at(i),
					0);
			initialVelocities.at(1)(i) = analyticVelocity2(supportPoints.at(i),
					0);
			;
		}
	}

	/**
	 * @short analytic solution of the Taylor-Green vortex, first component of velocity vector
	 */
	double analyticVelocity1(const dealii::Point<2>& x, double t) const {
		return sin(x(0)) * cos(x(1)) * exp(-2 * getViscosity() * t);
	}

	/**
	 * @short analytic solution of the Taylor-Green vortex, second component of velocity vector
	 */
	double analyticVelocity2(const dealii::Point<2>& x, double t) const {
		return -cos(x(0)) * sin(x(1)) * exp(-2 * getViscosity() * t);
	}

private:

	size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid() {
		//Creation of the principal domain
		boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(
#ifdef WITH_TRILINOS_MPI
				MPI_COMM_WORLD
#endif
				);
		dealii::GridGenerator::hyper_cube(*square, 0, 2 * Math::PI);

		// Assign boundary indicators to the faces of the "parent cell"
		Mesh<2>::active_cell_iterator cell = square->begin_active();
		cell->face(0)->set_all_boundary_ids(0);  // left
		cell->face(1)->set_all_boundary_ids(1);  // right
		cell->face(2)->set_all_boundary_ids(2);  // top
		cell->face(3)->set_all_boundary_ids(3);  // bottom


		return square;
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

		// Get the triangulation object (which belongs to the parent class).
		boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

		return boundaries;
	}

	virtual void refine(){
		getMesh()->refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<2>& ) {

	}
}
;

#endif /* TaylorGreenTest2D_H_ */
