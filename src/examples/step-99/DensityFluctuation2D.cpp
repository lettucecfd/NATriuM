/*
 * DensityFluctuation2D.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: kk
 */

#include "DensityFluctuation2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "problemdescription/PeriodicBoundary.h"

#include "utilities/Math.h"

namespace natrium {


DensityFluctuation2D::DensityFluctuation2D(double viscosity, size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(refinementLevel), viscosity, 1.0) {

	setBoundaries(makeBoundaries());

}

DensityFluctuation2D::~DensityFluctuation2D() {
}


void DensityFluctuation2D::applyInitialDensities(
		distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const {
	for (size_t i = 0; i < initialDensities.size(); i++) {
		if( (supportPoints.at(i)(0)<0.55)&&(supportPoints.at(i)(0)>0.45)
				&&(supportPoints.at(i)(1)<0.55)&&(supportPoints.at(i)(1)>0.45) ){
			initialDensities(i)=1.1 ;
		}
		else{
			initialDensities(i) = 1.0;
		}
	}
}

void DensityFluctuation2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const {
	assert(
			initialVelocities.at(0).size()
					== initialVelocities.at(1).size());
	for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
		initialVelocities.at(0)(i) = 0.0;
		initialVelocities.at(1)(i) = 0.0;
	}
}

shared_ptr<Triangulation<2> > DensityFluctuation2D::makeGrid(size_t refinementLevel) {
	shared_ptr<Triangulation<2> > square = make_shared<Triangulation<2> >();
	dealii::GridGenerator::hyper_cube(*square,0,1.);
	Triangulation<2>::active_cell_iterator cell = square->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // top
	cell->face(3)->set_all_boundary_indicators(3);  // bottom

	square->refine_global(refinementLevel);

	return square;
}

shared_ptr<BoundaryCollection<2> > DensityFluctuation2D::makeBoundaries() {
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<BoundaryCollection<2> >();
	boundaries->addBoundary(make_shared<PeriodicBoundary<2> >(0, 1, getTriangulation()));
	boundaries->addBoundary(make_shared<PeriodicBoundary<2> >(2, 3, getTriangulation()));
	shared_ptr<Triangulation<2> > tria_pointer = getTriangulation();
	return boundaries;
}
} /* namespace natrium */
