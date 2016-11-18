/**
 * @file PoiseuilleFlow3D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PoiseuilleFlow3D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/tensor.h"

#include "../boundaries/LinearFluxBoundaryRhoU.h"
#include "../boundaries/PeriodicBoundary.h"
#include "../problemdescription/ConstantExternalForce.h"
#include "../utilities/Math.h"

namespace natrium {

PoiseuilleFlow3D::PoiseuilleFlow3D(double viscosity, size_t refinementLevel,
		double u_bulk, double height, double width, double length, bool is_periodic) :
		Benchmark<3>(makeGrid(height, width, length), viscosity, height), m_uBulk(
				u_bulk), m_uMax(3. / 2. * u_bulk), m_refinementLevel(refinementLevel) {

	/// apply boundary values
	setBoundaries(makeBoundaries(is_periodic));
	// apply initial values / analytic solution
	setAnalyticU(boost::make_shared<AnalyticVelocity>(this));

	if (is_periodic) {
		// add external force
		double Fx = 8 * m_uMax * viscosity / (height * height);
		//pout << "F: " << Fx << endl;
		dealii::Tensor<1, 3> F;
		F[0] = Fx;
		setExternalForce(
				boost::make_shared<ConstantExternalForce<3> >(F));
	}

}

PoiseuilleFlow3D::~PoiseuilleFlow3D() {
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > PoiseuilleFlow3D::makeGrid(double height, double width,
		double length) {
	//Creation of the principal domain
	vector<unsigned int> repetitions(3);
	repetitions.at(0) = 1;
	repetitions.at(1) = 1;
	repetitions.at(2) = 1;

	boost::shared_ptr<Mesh<3> > mesh = boost::make_shared<Mesh<3> >(
	MPI_COMM_WORLD);

	bool colorize = true; 	// do not set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	dealii::GridGenerator::subdivided_hyper_rectangle(*mesh, repetitions,
			dealii::Point<3>(0.0, 0.0, 0.0),
			dealii::Point<3>(length, width, height), colorize);
	return mesh;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3> > PoiseuilleFlow3D::makeBoundaries(
		bool is_periodic) {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
			BoundaryCollection<3> >();
	dealii::Vector<double> zeroVector(3);

	if (is_periodic) {
//		boundaries->addBoundary(
//				boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
//		boundaries->addBoundary(
//				boost::make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
//		boundaries->addBoundary(
//				boost::make_shared<LinearBoundaryRhoU<3> >(4, zeroVector));
//		boundaries->addBoundary(
//				boost::make_shared<LinearBoundaryRhoU<3> >(5, zeroVector));
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh())); // TODO: damn, comment it!
		//cout << " > periodic: inlet/outlet" << endl;
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
		//cout << " > periodic: back/front" << endl;
		boundaries->addBoundary(
				boost::make_shared<LinearFluxBoundaryRhoU<3> >(4, zeroVector));
		//cout << " > no-slip: top" << endl;
		boundaries->addBoundary(
				boost::make_shared<LinearFluxBoundaryRhoU<3> >(5, zeroVector));
		//cout << " > no-slip: bottom" << endl;

	} else {
		natrium_errorexit("Periodic-free Poiseuille flow not implemented, yet.");
//		dealii::Vector<double> xVelocity(2);
//		xVelocity(0) = m_uMax;
//		boundaries->addBoundary(
//				boost::make_shared<LinearBoundaryRhoU<3> >(0,
//						boost::make_shared<dealii::ConstantFunction<3> >(1.0),
//						boost::make_shared<PoiseuilleFlow3D::AnalyticVelocity>(
//								this)));
//		boundaries->addBoundary(
//				boost::make_shared<NonlinearBoundaryZouHeRho<3> >(1,
//						boost::make_shared<dealii::ConstantFunction<3> >(-1.0),
//						1));
//		boundaries->addBoundary(
//				boost::make_shared<LinearBoundaryRhoU<3> >(2, zeroVector));
//		boundaries->addBoundary(
//				boost::make_shared<LinearBoundaryRhoU<3> >(3, zeroVector));
	}

	return boundaries;
}

double PoiseuilleFlow3D::AnalyticVelocity::value(const dealii::Point<3>& ,
		const unsigned int component) const {
	assert(component < 3);
	/*double h = m_flow->getCharacteristicLength();
	if (component == 0) {
		return (- 4 * m_flow->m_uMax *
				(x(2) - h) * x(2) / (h*h) );
	} else {*/
		return 0.0;
	//}
}

} /* namespace natrium */

