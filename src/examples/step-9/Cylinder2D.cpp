/**
 * @file Cylinder2D.cpp
 * @short 
 * @date 09.10.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Cylinder2D.h"

#include "deal.II/grid/grid_in.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

#include "natrium/problemdescription/PeriodicBoundary.h"
#include "natrium/problemdescription/MinLeeBoundary.h"

#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/Logging.h"

namespace natrium {


Cylinder2D::Cylinder2D(double viscosity, double inletVelocity) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_inletVelocity(inletVelocity) {
	setCharacteristicLength(1.0);

	/// apply boundary values
	setBoundaries(makeBoundaries(inletVelocity));
}

Cylinder2D::~Cylinder2D() {
}

boost::shared_ptr<Mesh<2> > Cylinder2D::makeGrid() {

	//Creation of the principal domain
	boost::shared_ptr<Mesh<2> > tria = boost::make_shared<Mesh<2> >();
	dealii::GridIn<2> gridin;
	gridin.attach_triangulation(*tria);
	// !!! When including unv meshes from Salome, the first two blocks have to be deleted manually (compare Original_Mesh_1.unv to Mesh_1.unv)
	// The mesh file then regularly begins with the lines "-1" and "2411"
	std::ifstream f("/home/kraemer/eclipse_workspace/NATriuM/src/examples/step-9/salome/Mesh_1.unv");
	gridin.read_unv(f);
	CFDSolverUtilities::mesh_info(*tria, "cylinder.eps");

	return tria;
}

boost::shared_ptr<BoundaryCollection<2> > Cylinder2D::makeBoundaries(
		double inletVelocity) {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = inletVelocity;

	boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(0, constantVelocity));
	boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(1, constantVelocity));
	boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(2, constantVelocity));
	boundaries->addBoundary(
			boost::make_shared<MinLeeBoundary<2> >(3, constantVelocity));
	boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(4, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}

/**
 * @short set initial densities
 * @param[out] initialDensities vector of densities; to be filled
 * @param[in] supportPoints the coordinates associated with each degree of freedom
 */
void  Cylinder2D::applyInitialDensities(distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const{
	for (size_t i = 0; i < initialDensities.size(); i++){
		initialDensities(i) = 1.0;
	}
}

/**
 * @short set initial velocities
 * @param[out] initialVelocities vector of velocities; to be filled
 * @param[in] supportPoints the coordinates associated with each degree of freedom
 */
 void  Cylinder2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const{
	 assert ( initialVelocities.size() == 2);
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++){
			initialVelocities.at(0)(i) = getCharacteristicVelocity();
			initialVelocities.at(1)(i) = 0.0;
		}
 }


} /* namespace natrium */
