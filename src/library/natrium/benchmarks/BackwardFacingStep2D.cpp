/*
 * BackwardFacingStep2D.cpp
 *
 *  Created on: Jan 05, 2016
 *      Author: kraemer
 */

#include "BackwardFacingStep2D.h"

#include <math.h>

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/grid_tools.h"

#include "../boundaries/LinearFluxBoundaryRhoU.h"
#include "../utilities/BasicNames.h"

namespace natrium {

BackwardFacingStep2D::BackwardFacingStep2D(double viscosity,
		double inflow_velocity, size_t refinement_level, double L_domain,
		double L_step, double h_domain, double h_step) :
		ProblemDescription<2>(makeGrid(), viscosity, h_domain), m_inflowVelocity(
				inflow_velocity), m_LDomain(L_domain), m_LStep(L_step), m_HDomain(
				h_domain), m_HStep(h_step), m_refinementLevel(refinement_level) {
	setBoundaries(makeBoundaries(inflow_velocity));
	this->setInitialRho(boost::make_shared<InitialVelocity>(this));

}

BackwardFacingStep2D::~BackwardFacingStep2D() {
}

boost::shared_ptr<Mesh<2> > BackwardFacingStep2D::makeGrid() {

	size_t nx_in = 2;
	size_t nx_out = 18;
	size_t ny_top = 1;
	size_t ny_bottom = 2;

	boost::shared_ptr<Mesh<2> > rect_lu = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
	boost::shared_ptr<Mesh<2> > rect_ru = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
	boost::shared_ptr<Mesh<2> > rect_rb = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);

	boost::shared_ptr<Mesh<2> > bfs = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
	std::vector<unsigned int> repetitions(2);
	bool colorize = false; 	// do not set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	repetitions.at(0) = nx_in;
	repetitions.at(1) = ny_top;
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect_lu, repetitions,
			dealii::Point<2>(0.0, m_HStep),
			dealii::Point<2>(m_LStep, m_HDomain), colorize);
	repetitions.at(0) = nx_out;
	repetitions.at(1) = ny_top;
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect_ru, repetitions,
			dealii::Point<2>(m_LStep, m_HStep),
			dealii::Point<2>(m_LDomain, m_HDomain), colorize);
	repetitions.at(0) = nx_out;
	repetitions.at(1) = ny_bottom;
	dealii::GridGenerator::subdivided_hyper_rectangle(*rect_lu, repetitions,
			dealii::Point<2>(m_LStep, 0.0),
			dealii::Point<2>(m_LDomain, m_HStep), colorize);

	// merge triangulations
	boost::shared_ptr<Mesh<2> > tmp = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
	dealii::GridGenerator::merge_triangulations(*rect_lu, *rect_ru, *tmp);
	dealii::GridGenerator::merge_triangulations(*tmp, *rect_rb, *bfs);

	// set boundary ids
	typename dealii::Triangulation<2>::active_cell_iterator cell =
			bfs->begin_active(), endc = bfs->end();
	for (; cell != endc; ++cell)
		for (unsigned int f = 0; f < dealii::GeometryInfo<2>::faces_per_cell;
				++f) {
			const dealii::Point<2> face_center = cell->face(f)->center();
			if (cell->face(f)->at_boundary()) {
				if (((fabs(face_center[1]) < 1e-6)
						or (fabs(face_center[1] - m_HDomain) < 1e-6))
						or (fabs(face_center[0] - m_LStep) < 1e-6)) {
					// solid wall
					cell->face(f)->set_boundary_id(0);
				} else if (fabs(face_center[0]) < 1e-6) {
					// inflow
					cell->face(f)->set_boundary_id(1);
				} else if (fabs(face_center[0] - m_LDomain) < 1e-6) {
					// outflow
					cell->face(f)->set_boundary_id(2);
				}
			}
		}

	return bfs;
}

boost::shared_ptr<BoundaryCollection<2> > BackwardFacingStep2D::makeBoundaries(
		double bottomVelocity) {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = bottomVelocity;
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<2> >(0, zeroVelocity));
	boundaries->addBoundary(
			boost::make_shared<LinearFluxBoundaryRhoU<2> >(1,
					boost::make_shared<dealii::ConstantFunction<2> >(0.0),
					boost::make_shared<BackwardFacingStep2D::InflowVelocity>(
							m_HStep, m_HDomain, m_inflowVelocity)));
	//TODO add pressure boundary
	/*boundaries->addBoundary(
			boost::make_shared<NonlinearBoundaryZouHeRho<2> >(2,
					boost::make_shared<dealii::ConstantFunction<2> >(0.0), 1));*/

	return boundaries;
}

void BackwardFacingStep2D::refine(Mesh<2>& mesh){
	// Refine grid
	mesh.refine_global(m_refinementLevel);
}
void BackwardFacingStep2D::transform(Mesh<2>& mesh){
	// transform grid
	//dealii::GridTools::transform(
	//		UnstructuredGridFunc(averageHeight, amplitude, L), mesh);
	std::ofstream out("grid-2.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(mesh, out);
}


double BackwardFacingStep2D::InflowVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component < 2);
	// Poiseuille profile
	if (component == 0)
		return (-m_uAverage * 1.5 * 4 * (x[1] - m_hDomain) * (x[1] - m_hStep)
				/ std::pow(m_hDomain - m_hStep, 2));
	else
		return 0.0;
}

double BackwardFacingStep2D::InitialVelocity::value(const dealii::Point<2>&,
		const unsigned int component) const {
	assert(component < 2);
#ifdef FASTER_CONVERGENCE_INIT
	if (component == 0) {
		double upper = m_flow->m_height + flow->m_ampl * std::sin(8 * std::atan(1) * x(0) );
		return flow->m_bottomVelocity * pow( 1 - x(1)/upper,2);
	}
#else
	return 0.0;
#endif
}

} /* namespace natrium */
