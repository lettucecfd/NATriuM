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

#include "../boundaries/PeriodicBoundary.h"
#include "../boundaries/VelocityNeqBounceBack.h"

#include "../utilities/Math.h"

namespace natrium {

SodShockTube::SodShockTube(double viscosity, size_t refinement_level, double u0,
		double kappa, double perturbation, double trafo_x, double trafo_y) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_u0(u0), m_kappa(
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

	cout << "Density initialized" << endl;
			if (x(0) <= 8.0) {
				return 4.0;
			}
		 else {
			return 1.0; }



}

double SodShockTube::InitialTemperature::value(const dealii::Point<2>& x, const unsigned int component) const {
	cout << "Temperature initialized" << endl;

			if (x(0) <= 8.0) {
				return 1.0;
			}
		 else {
			return 1.0; }



}



/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > SodShockTube::makeGrid() {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >();
#endif
	const dealii::Point<2> left = {0.0,0.0};
	const dealii::Point<2> right = {100.0,1.0};
	const std::vector <unsigned int>& reps = {100,1};


	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, reps, left, right, true);
	//dealii::GridGenerator::hyper_cube(*rect, 0, 1);
	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = rect->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom

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
