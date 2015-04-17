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

TaylorGreenVortex3D::TaylorGreenVortex3D(double viscosity, size_t refinementLevel) :
    Benchmark<3>(makeGrid(refinementLevel), viscosity, 1) {

  /// apply boundary values
  setBoundaries(makeBoundaries());

}

TaylorGreenVortex3D::~TaylorGreenVortex3D() {
}


void TaylorGreenVortex3D::getAnalyticVelocity(const dealii::Point<3>& x, double t, dealii::Point<3>& velocity) const {
	  velocity(0) = sin(x(0))*cos(x(1))*cos(x(2))*exp(-2*getViscosity()*t);
	  velocity(1) = -cos(x(0))*sin(x(1))*cos(x(2))*exp(-2*getViscosity()*t);
	  velocity(2) = 0;

}

/**
 * @short create triangulation for TaylorGreen Vortex flow
 * @return shared pointer to a triangulation instance
 */
shared_ptr<Triangulation<3> > TaylorGreenVortex3D::makeGrid(size_t refinementLevel) {
  //Creation of the principal domain
  shared_ptr<Triangulation<3> > square =
      make_shared<Triangulation<3> >();
  dealii::GridGenerator::hyper_cube(*square, 0, 8*atan(1));

  // Assign boundary indicators to the faces of the "parent cell"
  Triangulation<3>::active_cell_iterator cell =
      square->begin_active();
  cell->face(0)->set_all_boundary_indicators(0);  //
  cell->face(1)->set_all_boundary_indicators(1);  //
  cell->face(2)->set_all_boundary_indicators(2);  //
  cell->face(3)->set_all_boundary_indicators(3);  //
  cell->face(4)->set_all_boundary_indicators(4);  //
  cell->face(5)->set_all_boundary_indicators(5);  //

  // Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
  square->refine_global(refinementLevel);

  return square;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
shared_ptr<BoundaryCollection<3> > TaylorGreenVortex3D::makeBoundaries() {

  // make boundary description
  shared_ptr<BoundaryCollection<3> > boundaries = make_shared<
      BoundaryCollection<3> >();
  boundaries->addBoundary(
      make_shared<PeriodicBoundary<3> >(0, 1, getTriangulation()));
  boundaries->addBoundary(
      make_shared<PeriodicBoundary<3> >(2, 3, getTriangulation()));
  boundaries->addBoundary(
      make_shared<PeriodicBoundary<3> >(4, 5, getTriangulation()));

  // Get the triangulation object (which belongs to the parent class).
  shared_ptr<Triangulation<3> > tria_pointer = getTriangulation();

  return boundaries;
}
} /* namespace natrium */





