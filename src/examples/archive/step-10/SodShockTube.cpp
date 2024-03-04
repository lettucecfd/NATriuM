/**
 * @file SodShockTube.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SodShockTube.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/VelocityNeqBounceBack.h"

#include "natrium/utilities/Math.h"

namespace natrium {

SodShockTube::SodShockTube(int length, double viscosity, size_t refinement_level, double u0,
		double kappa, double perturbation, double trafo_x, double trafo_y) :
		  ProblemDescription<2>(makeGrid(length), viscosity, 1.0), m_length(length), m_u0(u0), m_kappa(
				kappa), m_refinementLevel(refinement_level), m_perturbation(perturbation),
       	              		m_trafoX(trafo_x), m_trafoY(trafo_y)	{
	assert(trafo_x >=0);
	assert(trafo_x < 1);
	assert(trafo_y >=0);
	assert(trafo_y < 1);
	/// apply boundary values
	setBoundaries(makeBoundaries());
	// apply initial and analytical solution
	//this->setInitialU(boost::make_shared<InitialVelocity>(this));
	this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

SodShockTube::~SodShockTube() {
}

double SodShockTube::InitialVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	return 0.0;

}

double SodShockTube::InitialDensity::value(const dealii::Point<2>& x, const unsigned int component) const {

    double middle = static_cast<double>(m_flow->m_length) / 2.0;
    double steepness = 10000.0;
    auto[tanh_c, tanh_d] = m_flow->calcOffsets(8.0, 1.0);

    double return_value = (tanh(steepness * (x[0] - middle)) + tanh_c) * tanh_d;

    return return_value;

}

double SodShockTube::InitialTemperature::value(const dealii::Point<2>& x, const unsigned int component) const {

    double middle = static_cast<double>(m_flow->m_length) / 2.0;
    double steepness = 10000.0;
    auto[tanh_c, tanh_d] = m_flow->calcOffsets(1.25, 1.0);

    double return_value = (tanh(steepness * (x[0] - middle)) + tanh_c) * tanh_d;

    return return_value;



}



/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > SodShockTube::makeGrid(int length) {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >();
#endif
	const dealii::Point<2> left = {0.0,0.0};
	const dealii::Point<2> right = {static_cast<double>(length),1.0};
	const std::vector <unsigned int>& reps = {static_cast<unsigned int>(length),1};


	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, reps, left, right, true);
	//dealii::GridGenerator::hyper_cube(*rect, 0, 1);
	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = rect->begin_active();
	//cell->face(0)->set_boundary_id(0);  // left
	//cell->face(1)->set_boundary_id(1);  // right
	//cell->face(2)->set_boundary_id(2);  // top
	//cell->face(3)->set_boundary_id(3);  // bottom

	return rect;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > SodShockTube::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
	numeric_vector zeroVelocity(2);

	boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(0, zeroVelocity));
	boundaries->addBoundary(
			boost::make_shared<VelocityNeqBounceBack<2> >(1, zeroVelocity));
	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));


	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
