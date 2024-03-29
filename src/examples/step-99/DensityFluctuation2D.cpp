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

#include "natrium/boundaries/PeriodicBoundary.h"

#include "natrium/utilities/Math.h"

namespace natrium {


DensityFluctuation2D::DensityFluctuation2D(double viscosity, size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(refinementLevel), viscosity, 1.0), m_refinementLevel(refinementLevel) {

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

boost::shared_ptr<Mesh<2> > DensityFluctuation2D::makeGrid() {
	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >();
	dealii::GridGenerator::hyper_cube(*square,0,1.);
	Mesh<2>::active_cell_iterator cell = square->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom

	return square;
}

boost::shared_ptr<BoundaryCollection<2> > DensityFluctuation2D::makeBoundaries() {
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<BoundaryCollection<2> >();
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();
	return boundaries;
}
} /* namespace natrium */
