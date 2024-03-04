/*
 * TaylorGreenVortex3D.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: bajat
 */

#include "TaylorGreenVortex3D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/grid_out.h"

#include "../boundaries/PeriodicBoundary.h"

namespace natrium {

TaylorGreenVortex3D::TaylorGreenVortex3D(double viscosity, size_t refinementLevel, double cs,
                                         bool init_rho_analytically, size_t repetitions, double densityNumerator,
                                         bool compressible_case) :
        ProblemDescription<3>(makeGrid(repetitions), viscosity, 1), m_cs(cs), m_densityNumerator(densityNumerator),
        m_analyticInit(init_rho_analytically), m_compressibleCase(compressible_case), m_refinementLevel(refinementLevel)
{
	/// apply boundary values
	setBoundaries(makeBoundaries());
	// apply analytic solution
	this->setInitialU(boost::make_shared<InitialVelocity>(this));
	this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

TaylorGreenVortex3D::~TaylorGreenVortex3D() {
}

double TaylorGreenVortex3D::InitialVelocity::value(const dealii::Point<3>& x,
		const unsigned int component) const {
	assert(component < 3);
	if (component == 0) {
		return sin(x(0)) * cos(x(1)) * cos(x(2));
	} else if (component == 1) {
		return -cos(x(0)) * sin(x(1)) * cos(x(2));
	} else {
		return 0;
	}
}

double TaylorGreenVortex3D::InitialDensity::value(const dealii::Point<3>& x,
		const unsigned int component) const {
	assert(component == 0);
	if (m_flow->m_analyticInit) {
		double p = 1.0 / 16. * (cos(2 * x(0)) + cos(2 * x(1))) * (cos (2* x(2)) + 2);
		if(m_flow->m_compressibleCase){
            return 1.0 + p * m_flow->m_densityNumerator;
		}
        return 1.0 + p / (m_flow->m_cs * m_flow->m_cs);
    } else {
		return 1.0;
	}
}

    double TaylorGreenVortex3D::InitialTemperature::value(const dealii::Point<3>& x, const unsigned int component) const {
        (void)x;
        assert(component == 0);
        if (m_flow->m_analyticInit) {
            return 1.0;
        } else {
            return 1.0;
        }
    }

/**
 * @short create triangulation for TaylorGreen Vortex flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > TaylorGreenVortex3D::makeGrid(size_t repetitions) {
	//Creation of the principal domain

	boost::shared_ptr<Mesh<3> > cube =
			boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);

	dealii::Point<3> corner1(0,0,0);
	dealii::Point<3> corner2(8 * atan(1),8 * atan(1),8 * atan(1));
	std::vector<unsigned int> rep;
	rep.push_back(repetitions);
	rep.push_back(repetitions);
	rep.push_back(repetitions);
	dealii::GridGenerator::subdivided_hyper_rectangle(*cube, rep, corner1, corner2, true);

	// Assign boundary indicators to the faces of the "parent cell"
	/*Mesh<3>::active_cell_iterator cell = cube->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  //
	cell->face(1)->set_all_boundary_ids(1);  //
	cell->face(2)->set_all_boundary_ids(2);  //
	cell->face(3)->set_all_boundary_ids(3);  //
	cell->face(4)->set_all_boundary_ids(4);  //
	cell->face(5)->set_all_boundary_ids(5);  //
*/
	return cube;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3> > TaylorGreenVortex3D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
			BoundaryCollection<3> >();
	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<3> >(4, 5, 2, getMesh()));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */

