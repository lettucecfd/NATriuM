/*
 * Droplet2D.cpp
 *
 *  Created on: Dec 7, 2015
 *      Author: kraemer
 */

#include "Droplet2D.h"

#include "math.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/grid_tools.h"

#include "../problemdescription/LinearBoundaryRhoU.h"

namespace natrium {

Droplet2D::Droplet2D(double viscosity, size_t refinementLevel, double length, double height,
		double rho_l, double rho_g, double W, double R0) :
		ProblemDescription<2>(makeGrid(length,height), viscosity, height), m_length(length),
		m_height(height), m_rhoL(rho_l), m_rhoG(rho_g), m_W(W), m_R0(R0) {
//	assert(rho_l > rho_g);
//	assert(W > 0);
//	assert(W < 0.25 * L);
//	assert(R0 > 0);
//	assert(2 * R0 + 2 * W < L);

	setBoundaries(makeBoundaries());
	setInitialRho(boost::make_shared<InitialDensity>(this));
	setInitialU(boost::make_shared<InitialVelocity>(this));

	// refine grid
	boost::shared_ptr<Mesh<2> > rect = getMesh();
	rect->refine_global(refinementLevel);

}

Droplet2D::~Droplet2D() {
}

boost::shared_ptr<Mesh<2> > Droplet2D::makeGrid(double length, double height) {
	//Creation of the principal domain

//	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
//	dealii::GridGenerator::hyper_cube(*square, 0, L);

	// Assign boundary indicators to the faces of the "parent cell"
//	Mesh<2>::active_cell_iterator cell = square->begin_active();
//	cell->face(0)->set_all_boundary_ids(0);  // left
//	cell->face(1)->set_all_boundary_ids(1);  // right
//	cell->face(2)->set_all_boundary_ids(2);  // top
//	cell->face(3)->set_all_boundary_ids(3);  // bottom

	//	return square;

	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_rectangle(*rect, dealii::Point<2>(0, 0), dealii::Point<2>(length, height), false);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = rect->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // bottom
	cell->face(3)->set_all_boundary_ids(3);  // top

	return rect ;
}

boost::shared_ptr<BoundaryCollection<2> > Droplet2D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();
/*	numeric_vector zeroVelocity(2);
	boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(0, zeroVelocity));
	boundaries->addBoundary(
			boost::make_shared<MinLeeBoundary<2> >(1, zeroVelocity));
	boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(2, zeroVelocity));
	boundaries->addBoundary(
			boost::make_shared<MinLeeBoundary<2> >(3, zeroVelocity));
*/
	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
//	boundaries->addBoundary(
//			boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));

	dealii::Vector<double> zeroVelocity(2);
	boundaries->addBoundary(boost::make_shared<LinearBoundaryRhoU<2> >(2, zeroVelocity)) ;
	boundaries->addBoundary(boost::make_shared<LinearBoundaryRhoU<2> >(3, zeroVelocity)) ;

//	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh() ;

	return boundaries;
}

double Droplet2D::InitialDensity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component == 0);
	double rho_l = m_flow->getRhoL();
	double rho_g = m_flow->getRhoG();
	double W = m_flow->getW();
	double R0 = m_flow->getR0();
	double length = m_flow->getLength() ;
	double height = m_flow->getHeight() ;

	double fluctuation = rho_l*(1e-4) ;

/*	if ( (x(0)>length/4) && (x(0)<length/2) && (x(1)>height/4) && (x(1)<height/2) ){
//		return rho_l + fluctuation ;
		return rho_l ;
	}
	else {
		return rho_l ;
	}
*/
	double dist = std::sqrt(
			std::pow(x(0) - 0.5*length, 2) + std::pow(x(1) - 0.5*height, 2));
	double arg_tanh = (dist - R0) / W ;

	return 0.5 * (rho_l + rho_g)
			- 0.5 * (rho_l - rho_g) * std::sinh(arg_tanh) / std::cosh(arg_tanh) ;


}

double Droplet2D::InitialVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component < 2);
	// Initialize with poiseuille velocity profile

	if (component == 0) {
//		return 4. * 0.0001 * (x(1)-(x(1)*x(1)));
		return 0. ;
	} else {
		return 0.0;
	}
}

} /* namespace natrium */
