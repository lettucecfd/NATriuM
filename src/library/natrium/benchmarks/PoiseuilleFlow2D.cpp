/**
 * @file PoiseuilleFlow2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PoiseuilleFlow2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "../problemdescription/DirichletBoundaryRhoU.h"

#include "../utilities/Math.h"

namespace natrium {

PoiseuilleFlow2D::PoiseuilleFlow2D(double viscosity, size_t refinementLevel,
		double u_bulk, double height, double length, bool is_periodic) :
		Benchmark<2>(makeGrid(height, length),
				viscosity, height), m_uBulk(u_bulk), m_uMax(3. / 2. * u_bulk) {

	/// apply boundary values
	setBoundaries(makeBoundaries(is_periodic));
	// apply initial values / analytic solution
	setAnalyticU(boost::make_shared<AnalyticVelocity>(this));

	// refine global
	getMesh()->refine_global(refinementLevel);
}

PoiseuilleFlow2D::~PoiseuilleFlow2D() {
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > PoiseuilleFlow2D::makeGrid(double height, double length) {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >();
#endif
	dealii::GridGenerator::hyper_rectangle(*rect, dealii::Point<2>(0, -height),
			dealii::Point<2>(length, height), false);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = rect->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // bottom
	cell->face(3)->set_all_boundary_ids(3);  // top


	return rect;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > PoiseuilleFlow2D::makeBoundaries(
		bool is_periodic) {

	if (is_periodic){

	} else {

	}
	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	dealii::Vector<double> zeroVector(2);
	dealii::Vector<double> xVelocity(2);
	xVelocity(0) = 0.1 / sqrt(3);
	boundaries->addBoundary(boost::make_shared<DirichletBoundaryRhoU<2> >(0, zeroVector));
	boundaries->addBoundary(boost::make_shared<DirichletBoundaryRhoU<2> >(1, zeroVector));
	boundaries->addBoundary(boost::make_shared<DirichletBoundaryRhoU<2> >(2, zeroVector));
	boundaries->addBoundary(boost::make_shared<DirichletBoundaryRhoU<2> >(3, xVelocity));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}

double PoiseuilleFlow2D::AnalyticVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert (component < 2);
	if (component == 0) {
		return m_flow->m_uMax
				* (1 - pow(x(1) / m_flow->getCharacteristicLength(), 2));
	} else {
		return 0.0;
	}
}

} /* namespace natrium */

