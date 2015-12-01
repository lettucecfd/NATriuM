/*
 * TaylorGreenVortex3D.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: bajat
 */

#include "TaylorGreenVortex3D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "../problemdescription/PeriodicBoundary.h"

namespace natrium {

TaylorGreenVortex3D::TaylorGreenVortex3D(double viscosity,
		size_t refinementLevel) :
		Benchmark<3>(makeGrid(), viscosity, 1) {

	/// apply boundary values
	setBoundaries(makeBoundaries());
	// apply analytic solution
	this->setAnalyticU(boost::make_shared<AnalyticVelocity>(this));

	// Refine grid
	getMesh()->refine_global(refinementLevel);

}

TaylorGreenVortex3D::~TaylorGreenVortex3D() {
}

double TaylorGreenVortex3D::AnalyticVelocity::value(const dealii::Point<3>& x,
		const unsigned int component) const {
	assert(component < 3);
	if (component == 0) {
		return sin(x(0)) * cos(x(1)) * cos(x(2))
				* exp(-2 * m_flow->getViscosity() * this->get_time());
	} else if (component == 1) {
		return -cos(x(0)) * sin(x(1)) * cos(x(2))
				* exp(-2 * m_flow->getViscosity() * this->get_time());
	} else {
		return 0;
	}
}

/**
 * @short create triangulation for TaylorGreen Vortex flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > TaylorGreenVortex3D::makeGrid() {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<3> > square =
	boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<3> > square = boost::make_shared<Mesh<3> >();
#endif
	dealii::GridGenerator::hyper_cube(*square, 0, 8 * atan(1));

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<3>::active_cell_iterator cell = square->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  //
	cell->face(1)->set_all_boundary_ids(1);  //
	cell->face(2)->set_all_boundary_ids(2);  //
	cell->face(3)->set_all_boundary_ids(3);  //
	cell->face(4)->set_all_boundary_ids(4);  //
	cell->face(5)->set_all_boundary_ids(5);  //

	return square;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3> > TaylorGreenVortex3D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
			BoundaryCollection<3> >();
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(4, 5, 2, getMesh()));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */

