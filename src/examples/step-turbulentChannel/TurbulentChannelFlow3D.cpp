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
		double U_in, double height, double length, double width, bool is_periodic ) :
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
						SHIFTING_VELOCITY));
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

	// Synthetic turbulence generation due to Davidson et al.
	// [1] Using isotropic...
	// [2]
	/*
	charLength -> h
	visc -> viscosity
	*/
	// turbulent kinetic energy to be set properly
	//

	assert(component < 3);

	//------------------------------------------------------------------------------------------------------------------------------------------------------
	// Selectable variables
	// TODO: These variables has to be made selectable from a property file or similar

	double 	h = m_flow->getCharacteristicLength();				// characteristic flow length scale, here: full channel height
	double	qm 					= 4.0;							// turbulent kinetic energy
	double	delta 				= 0.5*h;        				// inlet boundary layer thickness. Only if the boundary layer
																// at inlet is not fully developed
	//double 	uTau 				= 1/25.0;						// inlet shear velocity, suTau = Urms [Ref2]

	// blending function parameters
	double	blendDist 			= 0.1*h;						// distance over which fBlend goes from 0 to 1
	double	freeStreamTurb 		= 0.1;							// parameter does not let fBlend drop below the prescribed value

	int 	nmodes 				= 150;							// number of Fourier modes
	double	wew1fct				= 2;							// ratio of ke and kmin (in wavenumber)

	// Constants
	double 	amp 				= 1.452762113;					// alpha in [Ref2]

	// Calculated
	double	sli 				= 0.1*delta;    				// length scale
	double	up 					= sqrt(2*qm/3);					// turbulent velocity scale
	double	epsm 				= pow(qm,1.5)/sli;				// dissipation rate



	if (component == 0) {
		// exponential law profile
		if (x(2) <= h/2)
			return ( m_flow->m_uBulk * std::pow(2*x(2)/h, 1./7.) );
		else
			return ( m_flow->m_uBulk * std::pow(2*(1-x(2)/h), 1./7.) );
		/*
		// parabolic profile
		return (- 4 * m_flow->m_uBulk * 1.5 *
				(x(2) - h) * x(2) / (h*h) );
		*/
	} else {
		return 0.0;
	}
}

} /* namespace natrium */

