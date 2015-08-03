/*
 * LubricationSine.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: kraemer
 */

#include "LubricationSine.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/grid_tools.h"

namespace natrium {

LubricationSine::LubricationSine(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L, double averageHeight,
		double amplitude, double cellAspectRatio,  double roughnessHeight, size_t roughnessLengthRatio) :
		ProblemDescription<2>(
				makeGrid(L, refinementLevel, averageHeight, amplitude, cellAspectRatio, roughnessHeight, roughnessLengthRatio),
				viscosity, averageHeight), m_bottomVelocity(bottomVelocity), m_height(averageHeight), m_ampl(amplitude), m_length(L) {
	setBoundaries(makeBoundaries(bottomVelocity));
}

LubricationSine::~LubricationSine() {
}

shared_ptr<Mesh<2> > LubricationSine::makeGrid(double L,
		size_t refinementLevel, double averageHeight, double amplitude, double cellAspectRatio,  double roughnessHeight, size_t roughnessLengthRatio) {

	//Creation of the principal domain
	shared_ptr<Mesh<2> > rect = make_shared<Mesh<2> >();
	dealii::Point<2> x1(0,0);
	dealii::Point<2> x2(L,averageHeight);
	std::vector<unsigned int> repetitions;
	repetitions.push_back( L / averageHeight / cellAspectRatio);
	repetitions.push_back( 1 );
	bool colorize = true; 	// set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, repetitions, x1, x2, colorize);

	// refine grid
	rect->refine_global(refinementLevel);

	// transform grid
	dealii::GridTools::transform(UnstructuredGridFunc(averageHeight, amplitude, L, roughnessHeight, roughnessLengthRatio),
			*rect);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(*rect, out);
	return rect;
}

shared_ptr<BoundaryCollection<2> > LubricationSine::makeBoundaries(
		double bottomVelocity) {

	// make boundary description
	shared_ptr<BoundaryCollection<2> > boundaries = make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = bottomVelocity;

	boundaries->addBoundary(
			make_shared<PeriodicBoundary<2> >(0, 1, getMesh()));
	boundaries->addBoundary(
			make_shared<MinLeeBoundary<2> >(2, constantVelocity));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(3, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}

void LubricationSine::applyInitialDensities(
		distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const {
	for (size_t i = 0; i < initialDensities.size(); i++) {
		initialDensities(i) = 1.0;
	}
}

void LubricationSine::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const {
	assert ( initialVelocities.size() == 2);
	for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
		double upper = m_height - m_ampl; //+  m_ampl * std::sin(8 * std::atan(1) * supportPoints.at(i)(0) / m_length() );
		initialVelocities.at(0)(i) = m_bottomVelocity * pow( 1 - supportPoints.at(i)(1)/upper,2);
		initialVelocities.at(1)(i) = 0.0;
	}
}

} /* namespace natrium */
