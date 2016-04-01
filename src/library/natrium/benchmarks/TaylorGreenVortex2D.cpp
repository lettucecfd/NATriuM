/**
 * @file TaylorGreenVortex2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TaylorGreenVortex2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "../problemdescription/PeriodicBoundary.h"

#include "../utilities/Math.h"

namespace natrium {

TaylorGreenVortex2D::TaylorGreenVortex2D(double viscosity,
		size_t refinementLevel, double cs, bool init_rho_analytically) :
		Benchmark<2>(makeGrid(), viscosity, 8 * atan(1)), m_cs(
				cs), m_analyticInit(init_rho_analytically), m_refinementLevel(refinementLevel) {

	/// apply boundary values
	setBoundaries(makeBoundaries());
	// apply initial and analytical solution
	this->setAnalyticU(boost::make_shared<AnalyticVelocity>(this));
	this->setAnalyticRho(boost::make_shared<AnalyticDensity>(this));

}

TaylorGreenVortex2D::~TaylorGreenVortex2D() {
}

double TaylorGreenVortex2D::AnalyticVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component < 2);
	if (component == 0) {
		return sin(x(0)) * cos(x(1))
				* exp(-2 * m_flow->getViscosity() * this->get_time());
	} else {
		return -cos(x(0)) * sin(x(1))
				* exp(-2 * m_flow->getViscosity() * this->get_time());
	}
}

double TaylorGreenVortex2D::AnalyticDensity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert (component == 0);
	if (m_flow->m_analyticInit) {
		double rho0 = 1;
		double p = rho0 / 4. * (cos(2 * x(0)) + cos(2 * x(1)))
				* exp(-4 * m_flow->getViscosity() * this->get_time());
		return rho0 + p / (m_flow->m_cs * m_flow->m_cs);
	} else {
		return 1.0;
	}
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > TaylorGreenVortex2D::makeGrid() {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >();
#endif
	dealii::GridGenerator::hyper_cube(*square, 0, 8 * atan(1));

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
boost::shared_ptr<BoundaryCollection<2> > TaylorGreenVortex2D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
