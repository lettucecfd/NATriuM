/**
 * @file DiamondObstacle2D.cpp
 * @short Flow around a diamond-shaped obstacle
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DiamondObstacle2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/grid_in.h"

#include "natrium/boundaries/LinearFluxBoundaryRhoU.h"

#include "natrium/utilities/Math.h"

namespace natrium {



DiamondObstacle2D::DiamondObstacle2D(double velocity, double viscosity,
		size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(refinementLevel), viscosity, 1.0), m_meanInflowVelocity(
				velocity), m_refinementLevel(refinementLevel) {

	/// apply boundary values
	setBoundaries(makeBoundaries());

}

DiamondObstacle2D::~DiamondObstacle2D() {
}


/**
 * @short Read the mesh given in deal.ii's step-35
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > DiamondObstacle2D::makeGrid(
		size_t refinementLevel) {
	//Read in grid
	//Taken from step-35 in deal.ii
	dealii::GridIn<2> grid_in;
	boost::shared_ptr<Mesh<2> > mesh = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
	grid_in.attach_triangulation(*mesh);
	{
		std::stringstream filename;
		filename << getenv("NATRIUM_DIR") << "/src/examples/step-grid-in/nsbench2.inp";
		std::ifstream file(filename.str().c_str());
		assert(file);
		grid_in.read_ucd(file);
	}
	 mesh->refine_global (refinementLevel);
	 mesh->get_boundary_ids();

	return mesh;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > DiamondObstacle2D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	dealii::Vector<double> zeroVector(2);
	boost::shared_ptr<dealii::Function<2> > boundary_density = boost::make_shared<
			dealii::ConstantFunction<2> > (1.0);
	boost::shared_ptr<dealii::Function<2> > boundary_velocity = boost::make_shared<
			InflowVelocity> (m_meanInflowVelocity);
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<2> >(1, zeroVector));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<2> >(2, boundary_density, boundary_velocity));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<2> >(3, boundary_density, boundary_velocity));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<2> >(4, zeroVector));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
