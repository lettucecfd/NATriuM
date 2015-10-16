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

#include "../utilities/BasicNames.h"

namespace natrium {

SinusoidalShear2D::SinusoidalShear2D(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L, double averageHeight,
		double amplitude, double cell_aspect_ratio) :
		ProblemDescription<2>(
				makeGrid(L, averageHeight, cell_aspect_ratio),
				viscosity, averageHeight), m_bottomVelocity(bottomVelocity), m_height(averageHeight), m_ampl(amplitude) {
	setBoundaries(makeBoundaries(bottomVelocity));
	this->setInitialRho(make_shared<InitialVelocity>(this));

	// refine grid
	shared_ptr<Mesh<2> > rect = getMesh();
	rect->refine_global(refinementLevel);

	// transform grid
	dealii::GridTools::transform(UnstructuredGridFunc(averageHeight, amplitude, L),
			*rect);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(*rect, out);

}

SinusoidalShear2D::~SinusoidalShear2D() {
}

shared_ptr<Mesh<2> > SinusoidalShear2D::makeGrid(double L, double averageHeight, double cell_aspect_ratio) {
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
			make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(
			make_shared<MinLeeBoundary<2> >(2, constantVelocity));
	boundaries->addBoundary(make_shared<MinLeeBoundary<2> >(3, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}

double SinusoidalShear2D::InitialVelocity::value(const dealii::Point<2>& x, const unsigned int component) const{
	assert (component < 2);
#ifdef FASTER_CONVERGENCE_INIT
	if (component == 0){
		double upper = m_flow->m_height +  flow->m_ampl * std::sin(8 * std::atan(1) * x(0) );
		return flow->m_bottomVelocity * pow( 1 - x(1)/upper,2);
	}
#else
		return 0.0;
#endif
}

} /* namespace natrium */
