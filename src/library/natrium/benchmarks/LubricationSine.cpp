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
#include "../problemdescription/LinearBoundaryRhoU.h"

namespace natrium {

LubricationSine::LubricationSine(double viscosity, double bottomVelocity,
		size_t refinementLevel, double L, double averageHeight,
		double amplitude, double cellAspectRatio, double roughnessHeight,
		size_t roughnessLengthRatio) :
		ProblemDescription<2>(makeGrid(L, averageHeight, cellAspectRatio),
				viscosity, averageHeight), m_bottomVelocity(bottomVelocity), m_height(
				averageHeight), m_ampl(amplitude), m_length(L), m_roughnessHeight(
				roughnessHeight), m_roughnessLengthRatio(roughnessLengthRatio), m_refinementLevel(
				refinementLevel) {
	setBoundaries(makeBoundaries(bottomVelocity));
	setInitialU(boost::make_shared<InitialVelocity>(this));
}

LubricationSine::~LubricationSine() {
}

boost::shared_ptr<Mesh<2> > LubricationSine::makeGrid(double L,
		double averageHeight, double cellAspectRatio) {

	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >();
#endif
	dealii::Point<2> x1(0, 0);
	dealii::Point<2> x2(L, averageHeight);
	std::vector<unsigned int> repetitions;
	repetitions.push_back(L / averageHeight / cellAspectRatio);
	repetitions.push_back(1);
	bool colorize = true; 	// set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, repetitions, x1,
			x2, colorize);

	return rect;
}

boost::shared_ptr<BoundaryCollection<2> > LubricationSine::makeBoundaries(
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
			boost::make_shared<LinearBoundaryRhoU<2> >(2, constantVelocity));
	boundaries->addBoundary(
			boost::make_shared<LinearBoundaryRhoU<2> >(3, zeroVelocity));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}

double LubricationSine::InitialVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component < 2);
	if (component == 0) {
		double upper = m_flow->m_height - m_flow->m_ampl; //+  m_ampl * std::sin(8 * std::atan(1) * supportPoints.at(i)(0) / m_length() );
		return m_flow->m_bottomVelocity * pow(1 - x(1) / upper, 2);
	} else {
		return 0.0;
	}
}

void LubricationSine::refine(Mesh<2>& mesh) {
	// refine grid
	mesh.refine_global(m_refinementLevel);
}

void LubricationSine::transform(Mesh<2>& mesh){
	// transform grid
	dealii::GridTools::transform(
			UnstructuredGridFunc(m_height, m_ampl, m_length,
					m_roughnessHeight, m_roughnessLengthRatio), mesh);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(mesh, out);
}

} /* namespace natrium */
