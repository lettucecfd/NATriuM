/**
 * @file ComplexWall1.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "ComplexWall1.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/utilities/Logging.h"

namespace natrium {

ComplexWall1::ComplexWall1(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L) :
		ProblemDescription<2>(makeGrid(L, refinementLevel), viscosity, L), m_bottomVelocity(
				bottomVelocity), m_refinementLevel(refinementLevel), m_L(L) {
	setCharacteristicLength(L);

	/// apply boundary values
	setBoundaries (makeBoundaries(bottomVelocity));}

ComplexWall1::~ComplexWall1() {
}

boost::shared_ptr<Mesh<2> > ComplexWall1::makeGrid(double L,
		size_t refinementLevel) {

	//Creation of the principal domain
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >();
	dealii::Point<2> x1(0,0);
	dealii::Point<2> x2(L,1);
	std::vector<unsigned int> repetitions;
	repetitions.push_back(10);
	repetitions.push_back( 1 );
	bool colorize = true; 	// set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, repetitions, x1, x2, colorize);
	return rect;

}

void ComplexWall1::refine(Mesh<2>& mesh) {
	// refine grid
	mesh.refine_global(m_refinementLevel);
}

void ComplexWall1::transform(Mesh<2>& mesh){
	// transform grid
	dealii::GridTools::transform(UnstructuredGridFunc(L),
			mesh);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(mesh, out);
}

boost::shared_ptr<BoundaryCollection<2> > ComplexWall1::makeBoundaries(
		double bottomVelocity) {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = bottomVelocity;

	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(2, constantVelocity));
	boundaries->addBoundary(boost::make_shared<VelocityNeqBounceBack<2> >(3, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}

/**
 * @short set initial densities
 * @param[out] initialDensities vector of densities; to be filled
 * @param[in] supportPoints the coordinates associated with each degree of freedom
 */
void  ComplexWall1::applyInitialDensities(distributed_vector& initialDensities,
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
 void  ComplexWall1::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const{
	 assert ( initialVelocities.size() == 2);
		for (size_t i = 0; i < initialVelocities.at(0).size(); i++){
			initialVelocities.at(0)(i) = 0.0;
			initialVelocities.at(1)(i) = 0.0;
		}
 }

} /* namespace natrium */
