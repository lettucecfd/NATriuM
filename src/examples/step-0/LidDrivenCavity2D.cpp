/**
 * @file LidDrivenCavity2D.cpp
 * @short Lid-driven cavity with three static walls and one moving wall
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "LidDrivenCavity2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/grid_tools.h"

#include "natrium/problemdescription/LinearBoundaryRhoU.h"

#include "natrium/utilities/Math.h"

namespace natrium {

LidDrivenCavity2D::LidDrivenCavity2D(double velocity, double viscosity,
		size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), topPlateVelocity(
				velocity), m_refinementLevel(refinementLevel) {

	/// apply boundary values
	setBoundaries(makeBoundaries());
}

LidDrivenCavity2D::~LidDrivenCavity2D() {
}

void LidDrivenCavity2D::refine(Mesh<2>& mesh){
	// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
	mesh.refine_global(m_refinementLevel);
}

void LidDrivenCavity2D::transform(Mesh<2>& mesh){
	dealii::GridTools::transform(UnstructuredGridFunc(0.2, 0.9, 0.15, 0.85), mesh);
}


/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > LidDrivenCavity2D::makeGrid() {
	//Creation of the principal domain
	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_rectangle(*square, dealii::Point<2>(0, 0),
			dealii::Point<2>(1, 1), false);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = square->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // bottom
	cell->face(3)->set_all_boundary_ids(3);  // top



	return square;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > LidDrivenCavity2D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	dealii::Vector<double> zeroVector(2);
	dealii::Vector<double> xVelocity(2);
	xVelocity(0) = topPlateVelocity;
	boundaries->addBoundary(boost::make_shared<LinearBoundaryRhoU<2> >(0, zeroVector));
	boundaries->addBoundary(boost::make_shared<LinearBoundaryRhoU<2> >(1, zeroVector));
	boundaries->addBoundary(boost::make_shared<LinearBoundaryRhoU<2> >(2, zeroVector));
	boundaries->addBoundary(boost::make_shared<LinearBoundaryRhoU<2> >(3, xVelocity));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
