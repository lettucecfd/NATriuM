/**
 * @file ShearLayer2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "ShearLayer2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "../problemdescription/PeriodicBoundary.h"

#include "../utilities/Math.h"

namespace natrium {

ShearLayer2D::ShearLayer2D(double viscosity, size_t refinement_level, double u0,
		double kappa) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_u0(u0), m_kappa(
				kappa), m_refinementLevel(refinement_level) {

	/// apply boundary values
	setBoundaries(makeBoundaries());
	// apply initial and analytical solution
	this->setInitialU(boost::make_shared<InitialVelocity>(this));

}

ShearLayer2D::~ShearLayer2D() {
}

double ShearLayer2D::InitialVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component < 2);
	if (component == 0) {
		if (x(1) <= 0.5) {
			return m_flow->m_u0 * tanh( m_flow->m_kappa * (x(1) - 0.25));
		} else {
			return  m_flow->m_u0 * tanh( m_flow->m_kappa * (0.75 - x(1)));
		}
	} else {
		double delta = 0.01;
		return delta *  m_flow->m_u0 * sin( 8 * atan(1) * (x(0) + 0.25));
	}
}


/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > ShearLayer2D::makeGrid() {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >();
#endif
	dealii::GridGenerator::hyper_cube(*square, 0, 1);

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
boost::shared_ptr<BoundaryCollection<2> > ShearLayer2D::makeBoundaries() {

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
} /* namespace natrium */
