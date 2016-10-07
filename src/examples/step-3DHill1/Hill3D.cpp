/**
 * @file Hill3D.cpp
 * @short Flow around a diamond-shaped obstacle
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/grid_in.h"

#include "natrium/utilities/Math.h"
#include "Hill3D.h"
#include "natrium/boundaries/LinearFluxBoundaryRhoU.h"

namespace natrium {



Hill3D::Hill3D(double velocity, double viscosity,
		size_t refinementLevel) :
		ProblemDescription<3>(makeGrid(refinementLevel), viscosity, 1.0), m_meanInflowVelocity(
				velocity), m_refinementLevel(refinementLevel) {

	/// apply boundary values
	setBoundaries(makeBoundaries());


}

Hill3D::~Hill3D() {
}


/**
 * @short Read the mesh given in deal.ii's step-35
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > Hill3D::makeGrid(
		size_t refinementLevel) {
	//Read in grid
	//Taken from step-35 in deal.ii
	dealii::GridIn<3> grid_in;
	boost::shared_ptr<Mesh<3> > mesh = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
	grid_in.attach_triangulation(*mesh);
	{
		// Gambit meshes can not be read -- has to be converted
		// e.g. to openfoam with OpenFoams fluent3DMeshToFoam utility
		// and then foamToTecplot360 (as an example)
		// However, it will have to be checked whether this conversion retains the boundary ids.
		std::stringstream filename;
		filename << getenv("NATRIUM_DIR") << "/src/examples/step-3DHill1/3DHill_testMeshFLUENT.msh";
		std::ifstream file(filename.str().c_str());
		assert(file);
		grid_in.read_unv(file);
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
boost::shared_ptr<BoundaryCollection<3> > Hill3D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
			BoundaryCollection<3> >();
	dealii::Vector<double> zeroVector(3);
	boost::shared_ptr<dealii::Function<3> > boundary_density = boost::make_shared<
			dealii::ConstantFunction<3> > (1.0);
	boost::shared_ptr<dealii::Function<3> > boundary_velocity = boost::make_shared<
			InflowVelocity> (m_meanInflowVelocity);
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<3> >(1, zeroVector));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<3> >(2, boundary_density, boundary_velocity));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<3> >(3, boundary_density, boundary_velocity));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<3> >(4, zeroVector));

	return boundaries;
}
} /* namespace natrium */
