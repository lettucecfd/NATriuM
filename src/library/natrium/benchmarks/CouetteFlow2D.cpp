/**
 * @file CouetteFlow2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CouetteFlow2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/geometry_info.h"

#include "../problemdescription/PeriodicBoundary.h"
#include "../problemdescription/LinearBoundaryRhoU.h"

#include "../utilities/Logging.h"
#include "../utilities/MPIGuard.h"

namespace natrium {

CouetteFlow2D::CouetteFlow2D(double viscosity, double topPlateVelocity,
		size_t refinementLevel, double L, double startTime, bool isUnstructured) :
		Benchmark<2>(makeGrid(L), viscosity, L), m_topPlateVelocity(
				topPlateVelocity), m_startTime(startTime), m_refinementLevel(refinementLevel), m_isUnstructured(isUnstructured) {
	setCharacteristicLength(L);

	//applyInitialValues
	this->setAnalyticU(boost::make_shared<AnalyticVelocity>(this));

	/// apply boundary values
	setBoundaries(makeBoundaries(topPlateVelocity));

}

CouetteFlow2D::~CouetteFlow2D() {
}

void CouetteFlow2D::refineAndTransform(){
	// refine grid
	boost::shared_ptr<Mesh<2> > unitSquare = getMesh();
	unitSquare->refine_global(m_refinementLevel);

	// transform grid
	if (m_isUnstructured) {
		dealii::GridTools::transform(UnstructuredGridFunc(), *unitSquare);
	}

	std::ofstream out("grid-couette.eps");
	dealii::GridOut grid_out;
	grid_out.write_eps(*unitSquare, out);
}


boost::shared_ptr<Mesh<2> > CouetteFlow2D::makeGrid(double L) {

	//Creation of the principal domain
	boost::shared_ptr<Mesh<2> > unitSquare = boost::make_shared<Mesh<2> >(
			MPI_COMM_WORLD);
	dealii::GridGenerator::hyper_cube(*unitSquare, 0, L);

	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = unitSquare->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // bottom
	cell->face(3)->set_all_boundary_ids(3);  // top

	return unitSquare;
}

boost::shared_ptr<BoundaryCollection<2> > CouetteFlow2D::makeBoundaries(
		double topPlateVelocity) {
	/*
	 // make boundary description
	 boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
	 BoundaryCollection<2> >();
	 numeric_vector zeroVelocity(2);
	 numeric_vector constantVelocity(2);
	 constantVelocity(0) = topPlateVelocity;

	 boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	 boundaries->addBoundary(boost::make_shared<MinLeeBoundary<2> >(2, zeroVelocity));
	 boundaries->addBoundary(
	 boost::make_shared<MinLeeBoundary<2> >(3, constantVelocity));

	 // Get the triangulation object (which belongs to the parent class).
	 boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();
	 */
	// make boundary description
	boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
			BoundaryCollection<2> >();

	numeric_vector zeroVelocity(2);
	numeric_vector constantVelocity(2);
	constantVelocity(0) = topPlateVelocity;

	boundaries->addBoundary(
			boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
	boundaries->addBoundary(
			boost::make_shared<LinearBoundaryRhoU<2> >(2, zeroVelocity));
	boundaries->addBoundary(
			boost::make_shared<LinearBoundaryRhoU<2> >(3, constantVelocity));


	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();
	return boundaries;
}

double CouetteFlow2D::AnalyticVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
	assert(component < 2);
	// the analytic solution is given by an asymptotic series
	if (component == 0) {
		double t = this->get_time();
		double U = m_benchmark->getCharacteristicVelocity();
		double L = m_benchmark->getCharacteristicLength();
		double t0 = m_benchmark->getStartTime();
		double nu = m_benchmark->getViscosity();

		t += t0;
		// the series converges veeeeeery slowly for t -> 0, thus assert t > epsilon
		// therefor the initial condition is set first:
		if (t < 0.00001) {
			if (x(1) < L - 0.00001) {
				return 0;
			} else {
				// upper border
				return U;
			}
		}

		double sum = U * x(1) / L;
		double lambda = 0.0;
		const double PI = atan(1) * 4;
		double increment = 1.0;
		double exp_expression = 0.0;

		for (size_t i = 1; i <= 10000; i++) {
			// calculate term in series
			lambda = i * PI / L;
			exp_expression = exp(-nu * lambda * lambda * t);
			increment = (i % 2 == 0 ? 1. : -1.) * 2 * U / (lambda * L)
					* exp_expression * sin(lambda * x(1));
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
	}
	return 0.0;
}

} /* namespace natrium */
