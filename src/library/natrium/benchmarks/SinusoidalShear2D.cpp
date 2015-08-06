/*
 * SinusoidalShear2D.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: kraemer
 */

#include "SinusoidalShear2D.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/grid_tools.h"

namespace natrium {

SinusoidalShear2D::SinusoidalShear2D(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L, double averageHeight,
		double amplitude, double cell_aspect_ratio) :
		ProblemDescription<2>(
				makeGrid(L, refinementLevel, averageHeight, amplitude, cell_aspect_ratio),
				viscosity, averageHeight), m_bottomVelocity(bottomVelocity), m_height(averageHeight), m_ampl(amplitude) {
	setBoundaries(makeBoundaries(bottomVelocity));
}

SinusoidalShear2D::~SinusoidalShear2D() {
}

shared_ptr<Mesh<2> > SinusoidalShear2D::makeGrid(double L,
		size_t refinementLevel, double averageHeight, double amplitude, double cell_aspect_ratio) {
#ifdef WITH_TRILINOS_MPI
	shared_ptr<Mesh<2> > rect = make_shared<Mesh<2> >(MPI_COMM_WORLD);
#else
	shared_ptr<Mesh<2> > rect = make_shared<Mesh<2> >();
#endif
	dealii::Point<2> x1(0,0);
	dealii::Point<2> x2(1,1);
	std::vector<unsigned int> repetitions;
	repetitions.push_back( 1./ cell_aspect_ratio);
	repetitions.push_back( 1 );
	bool colorize = true; 	// set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, repetitions, x1, x2, colorize);

	// refine grid
	rect->refine_global(refinementLevel);

	// transform grid
	dealii::GridTools::transform(UnstructuredGridFunc(averageHeight, amplitude, L),
			*rect);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(*rect, out);
	return rect;
}

shared_ptr<BoundaryCollection<2> > SinusoidalShear2D::makeBoundaries(
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

void SinusoidalShear2D::applyInitialDensities(
		distributed_vector& initialDensities,
		const vector<dealii::Point<2> >& supportPoints) const {
	for (size_t i = 0; i < initialDensities.size(); i++) {
		initialDensities(i) = 1.0;
	}
}

void SinusoidalShear2D::applyInitialVelocities(
		vector<distributed_vector>& initialVelocities,
		const vector<dealii::Point<2> >& supportPoints) const {
	assert ( initialVelocities.size() == 2);
	for (size_t i = 0; i < initialVelocities.at(0).size(); i++) {
#ifdef FASTER_CONVERGENCE_INIT
		double upper = m_height +  m_ampl * std::sin(8 * std::atan(1) * supportPoints.at(i)(0) );
		initialVelocities.at(0)(i) = m_bottomVelocity * pow( 1 - supportPoints.at(i)(1)/upper,2);
		initialVelocities.at(1)(i) = 0.0;
#else
		initialVelocities.at(0)(i) = 0.0;
		initialVelocities.at(1)(i) = 0.0;
#endif
	}
}

} /* namespace natrium */
