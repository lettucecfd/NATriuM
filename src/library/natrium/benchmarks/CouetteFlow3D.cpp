/**
 * @file CouetteFlow3D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CouetteFlow3D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

#include "../problemdescription/PeriodicBoundary.h"
#include "../problemdescription/MinLeeBoundary.h"

#include "../utilities/Logging.h"

namespace natrium {

CouetteFlow3D::CouetteFlow3D(double viscosity, double topPlateVelocity,
		size_t refinementLevel, double L, double startTime, bool isUnstructured) :
		Benchmark<3>(makeGrid(L, refinementLevel, isUnstructured), viscosity,
				L), m_topPlateVelocity(topPlateVelocity), m_startTime(startTime) {
	setCharacteristicLength(L);

	/// apply boundary values
	setBoundaries(makeBoundaries(topPlateVelocity));

	/// apply initial values
	setAnalyticU(make_shared<AnalyticVelocity>(this));
}

CouetteFlow3D::~CouetteFlow3D() {
}

shared_ptr<Mesh<3> > CouetteFlow3D::makeGrid(double L, size_t refinementLevel,
		bool isUnstructured) {

	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	shared_ptr<Mesh<3> > unitSquare = make_shared<Mesh<3> >(MPI_COMM_WORLD);
#else
	shared_ptr<Mesh<3> > unitSquare = make_shared<Mesh<3> >();
#endif
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, L);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<3>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_indicators(0);  // left
	cell->face(1)->set_all_boundary_indicators(1);  // right
	cell->face(2)->set_all_boundary_indicators(2);  // front
	cell->face(3)->set_all_boundary_indicators(3);  // back
	cell->face(4)->set_all_boundary_indicators(4);  // bottom
	cell->face(5)->set_all_boundary_indicators(5);  // top

	// refine grid
	unitSquare->refine_global(refinementLevel);

	// transform grid
	if (isUnstructured) {
		dealii::GridTools::transform(UnstructuredGridFunc(), *unitSquare);
	}
	return unitSquare;
}

shared_ptr<BoundaryCollection<3> > CouetteFlow3D::makeBoundaries(
		double topPlateVelocity) {

	// make boundary description
	shared_ptr<BoundaryCollection<3> > boundaries = make_shared<
			BoundaryCollection<3> >();
	numeric_vector zeroVelocity(3);
	numeric_vector constantVelocity(3);
	constantVelocity(0) = topPlateVelocity;

	boundaries->addBoundary(make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
	boundaries->addBoundary(make_shared<MinLeeBoundary<3> >(4, zeroVelocity));
	boundaries->addBoundary(
			make_shared<MinLeeBoundary<3> >(5, constantVelocity));

	// Get the triangulation object (which belongs to the parent class).
	shared_ptr<Mesh<3> > tria_pointer = getMesh();

	return boundaries;
}

double CouetteFlow3D::AnalyticVelocity::value(const dealii::Point<3>& x,
		const unsigned int component) const {
	assert(component < 3);
	if (component == 0) {
		// the analytic solution is given by an asymptotic series
		double t = this->get_time();
		double U = m_benchmark->getCharacteristicVelocity();
		double L = m_benchmark->getCharacteristicLength();
		double t0 = m_benchmark->getStartTime();
		double nu = m_benchmark->getViscosity();

		t += t0;
		// the series converges veeeeeery slowly for t -> 0, thus assert t > epsilon
		// therefor the initial condition is set first:
		if (t < 0.00001) {
			if (x(2) < L - 0.00001) {
				return 0;
			} else {
				// upper border
				return U;
			}
		}

		double sum = U * x(2) / L;
		double lambda = 0.0;
		const double PI = atan(1) * 4;
		double increment = 1.0;
		double exp_expression = 0.0;

		for (size_t i = 1; i <= 10000; i++) {
			// calculate term in series
			lambda = i * PI / L;
			exp_expression = exp(-nu * lambda * lambda * t);
			increment = (i % 2 == 0 ? 1. : -1.) * 2 * U / (lambda * L)
					* exp_expression * sin(lambda * x(2));
			// (i % 2 == 0 ? 1. : -1.) is a more efficient expression of (-1)^i
			sum += increment;
			// stop conditions: a) converged, b)
			if (exp_expression < 1e-20) {
				break;
			}
		}
		// assert convergence of the above asymptotic sum
		if (exp_expression >= 1e-12) {
			LOG(WARNING)
					<< "Warning: Analytic solution series did not converge."
					<< endl;
		}
		return sum;
	} else {
		return 0.0;
	}
}

} /* namespace natrium */
