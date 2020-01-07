/**
 * @file TurbulentChannelFlow3D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DecayingTurbulence2D.h"

#include <random>

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/tensor.h"
#include "deal.II/base/utilities.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include "natrium/problemdescription/ConstantExternalForce.h"
#include "natrium/utilities/Math.h"

namespace natrium {

DecayingTurbulence2D::DecayingTurbulence2D(double viscosity,
		size_t refinementLevel) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_refinementLevel(refinementLevel){

	/// apply boundary values
	setBoundaries (makeBoundaries());
	// apply initial values / analytic solution
	setInitialU	(boost::make_shared<InitialVelocity>(this));

}

DecayingTurbulence2D::~DecayingTurbulence2D() {
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > DecayingTurbulence2D::makeGrid() {

	boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(
			MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_cube(*square, 0, 1);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = square->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom

	return square;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > DecayingTurbulence2D::makeBoundaries() {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();

	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));

	return boundaries;
}

double DecayingTurbulence2D::InitialVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {

	assert(component < 2);
	double u = 0;
	double v = 0;
	double n, m, k;
	// spectral density function of Energy and Stream function
	double E, psi;
	double sine;
	// max wave number
	size_t max_k = 30;
	// spectral density function parameters
	double C0 = 1e-9;
	size_t A = 6;
	size_t B = 17;
	size_t k0 = 4;
	double PI = 4*std::atan(1);

	// make random divergence free fourier series
	for (size_t i = 0; i < 100; i++) {
		// set random seed, dependent on i
		srand(i+1000);
		// select random n
		//n = 1./ ( k0 * std::sqrt(2*PI) ) * std::exp(0.5 * std::pow( ((double) rand() / RAND_MAX) / k0, 2));
		n = dealii::Utilities::generate_normal_random_number(0, k0);
		//n = 2 + rand() % 3 ;
		// select random m
		//m = 1./ ( k0 * std::sqrt(2*PI) ) * std::exp(0.5 * std::pow( ((double) rand() / RAND_MAX) / k0, 2));
		//m = 5 + rand() % 30 ;
		m = dealii::Utilities::generate_normal_random_number(0, k0);
		k = std::sqrt(n*n + m*m);
		// sample random au_nm
		E = C0 * std::pow(k,A) / (1.0 + std::pow(k/k0,B+1) );
		pout << n << " " << m << " " << E << endl;
		psi = E / (k*k); // (spectral density of the stream function
		sine = std::sin(2*PI*n* x(0) + 2*PI*m* x(1) );
		// only sinus part
		u += (-psi * 2 * PI * m * sine);
		v += (psi * 2 * PI * n * sine);
	}
	pout << "u " << u << ", v " << v << endl;
	// return result;
	if (component == 0)
		return u;
	else
		return v;

}

} /* namespace natrium */

