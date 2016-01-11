/**
 * @file TurbulentChannelFlow3D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TurbulentChannelFlow3D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/base/tensor.h"

#include "natrium/problemdescription/LinearBoundaryRhoU.h"
#include "natrium/problemdescription/PeriodicBoundary.h"
#include "natrium/problemdescription/NonlinearBoundaryZouHeRho.h"
#include "natrium/problemdescription/ConstantExternalForce.h"
#include "natrium/utilities/Math.h"

namespace natrium {

TurbulentChannelFlow3D::TurbulentChannelFlow3D(double viscosity, size_t refinementLevel, double u_bulk,
		double height, double length, double width, bool is_periodic ) :
		ProblemDescription<3>(makeGrid(height, length, width), viscosity, height), m_uBulk(
				u_bulk) {

	/// apply boundary values
	setBoundaries(makeBoundaries(is_periodic));
	// apply initial values / analytic solution
	setInitialU(boost::make_shared<InitialVelocity>(this));

	if (is_periodic) {
		// add external force
		double Fx = 8 * 1.5 * m_uBulk * viscosity / (height * height);
		//pout << "F: " << Fx << endl;
		dealii::Tensor<1, 3> F;
		F[0] = Fx;
		setExternalForce(
				boost::make_shared<ConstantExternalForce<3> >(F,
						GUO));
	}

	// refine global
	getMesh()->refine_global(refinementLevel);
}

TurbulentChannelFlow3D::~TurbulentChannelFlow3D() {
}

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > TurbulentChannelFlow3D::makeGrid(double height,
		double length, double width) {
	//Creation of the principal domain

	boost::shared_ptr<Mesh<3> > mesh = boost::make_shared<Mesh<3> >(
	MPI_COMM_WORLD);

	std::vector<unsigned int> repetitions(3);
	bool colorize = true; 	// do not set boundary ids automatically to
							// 0:left; 1:right; 2:bottom; 3:top
	repetitions.at(0) = length / height;
	repetitions.at(1) = width / height;
	repetitions.at(2) = height;
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
boost::shared_ptr<BoundaryCollection<3> > TurbulentChannelFlow3D::makeBoundaries(
		bool is_periodic) {

	// make boundary description
	boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
			BoundaryCollection<3> >();
	dealii::Vector<double> zeroVector(3);

	if (is_periodic) {
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(4, zeroVector));
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(5, zeroVector));
	} else {
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(0,
						boost::make_shared<dealii::ConstantFunction<3> >(1.0),
						boost::make_shared<TurbulentChannelFlow3D::InitialVelocity>(
								this)));
		boundaries->addBoundary(
				boost::make_shared<NonlinearBoundaryZouHeRho<3> >(1,
						boost::make_shared<dealii::ConstantFunction<3> >(-1.0),
						1));
		boundaries->addBoundary(
						boost::make_shared<PeriodicBoundary<3> >(2, 3, 1, getMesh()));
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(4, zeroVector));
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(5, zeroVector));
	}

	return boundaries;
}

double TurbulentChannelFlow3D::InitialVelocity::value(const dealii::Point<3>& x,
		const unsigned int component) const {

	// TODO add synthetic turbulence
	assert(component < 3);
	double h = m_flow->getCharacteristicLength();
	if (component == 0) {
		return (- 4 * m_flow->m_uBulk * 1.5 *
				(x(2) - h) * x(2) / (h*h) );
	} else {
		return 0.0;
	}
}

} /* namespace natrium */

