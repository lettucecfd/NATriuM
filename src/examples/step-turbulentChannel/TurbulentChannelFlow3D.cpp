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

TurbulentChannelFlow3D::TurbulentChannelFlow3D(double viscosity, size_t refinementLevel,
		double u_bulk, double U_in, double height, double length, double width,
		double orderOfFiniteElement, bool is_periodic ) :
		ProblemDescription<3>(makeGrid(height, length, width), viscosity, height),
		m_refinementLevel(refinementLevel), m_uBulk(u_bulk), m_Uin(U_in),
		m_OFE(orderOfFiniteElement){

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
	//TODO: temporary variable
	//dealii::Tensor<1, 3> F = m_flow->getExternalForce();
	//cout << "External force" << F << endl;
	double ReTau = 7300;

	assert(component < 3);

	//------------------------------------------------------------------------------------------------------------------------------------------------------
	// Selectable variables
	// TODO: These variables have to be made selectable from a property file or similar

	double 	h = m_flow->getCharacteristicLength();				// characteristic flow length scale, here: full channel height
	double  visc = m_flow->getViscosity();						// kinematic viscosity
	double	qm 					= 4.0;							// turbulent kinetic energy
	double	delta 				= 5./32*h;	        			// inlet boundary layer thickness. Only if the boundary layer
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

	double 	OFE					= m_flow->getOrderOfFiniteElement();
	double	inletCellLength		= h/( (OFE + 1) * pow(2, m_flow->getRefinementLevel()) );
	double	wnrn				= 2*M_PI/inletCellLength;		// highest wave number // min cell length
	double	wnre				= 9*M_PI*amp/(55*sli);			// k_e (related to peak energy wave number)
	double	wnreta 				= pow((epsm/pow(visc, 3)), .25);
																// wavenumber used in the viscous expression (high wavenumbers)
																// 	in the von Karman spectrum
	double	wnr1 				= wnre/wew1fct;					// smallest wavenumber
	double	dkl 				= (wnrn-wnr1)/nmodes;			// wavenumber step

	// Create a seed from current time (must be negative).
	// 	A seed is a random input value for the algorithm generating
	// 	the random angles and is created once the code has started.
	int 	iseed	= -50;

	vector<double> 		fi(nmodes), psi(nmodes), alfa(nmodes), teta(nmodes);
	vector<int> 		iv(32);		// eq. to ntab in randf_2()
	int 				iy = 0;

	//TODO: show case details

	vector<double>		wnrf(nmodes), wnr(nmodes), dkn(nmodes);
	vector<double>		kxio(nmodes), kyio(nmodes), kzio(nmodes);
	vector<double>		sxio(nmodes), syio(nmodes), szio(nmodes);

	double				utrp, vtrp, wtrp;							// turbulent perturbation (fluctuation) in x, y and z direction
	double				kxi, kyi, kzi;
	double				sx, sy, sz, kx, ky, kz, rk;
	double				arg, tfunk, e, utn;

	double				minDist, fBlend;		 					// distance to the nearest wall

	//------------------------------------------------------------------------------------------------------------------------------------------------------
	// Wavenumber generation

	// Compute random angles
	m_flow->angles(nmodes, iseed, iy, iv, fi, psi, alfa, teta);

// DEBUG
//	for (int i = 0; i < 10; i++)
//		cout << fi[i] << endl;

	// Wavenumber at faces
	for (int m = 0; m <= nmodes; ++m)
	{
	wnrf[m] = wnr1+dkl*m;
	}

	// Wavenumber at cell centers
	for (int m = 0; m <= nmodes-1; ++m)
	{
	wnr[m] = 0.5*(wnrf[m+1]+wnrf[m]);
	dkn[m] = wnrf[m+1]-wnrf[m];
	}

	// Wavenumber vector from random angles
	for (int m = 0; m <= nmodes-1; ++m)
	{
	kxio[m]=sin(teta[m])*cos(fi[m]);
	kyio[m]=sin(teta[m])*sin(fi[m]);
	kzio[m]=cos(teta[m]);

	// sigma (s=sigma) from random angles. sigma is the unit direction which gives the direction
	// of the synthetic velocity vector (u, v, w)
	  sxio[m]=cos(fi[m])*cos(teta[m])*cos(alfa[m])-sin(fi[m])*sin(alfa[m]);
	  syio[m]=sin(fi[m])*cos(teta[m])*cos(alfa[m])+cos(fi[m])*sin(alfa[m]);
	  szio[m]=-sin(teta[m])*cos(alfa[m]);
	}

//	#makeFieldInc
//	//TODO: declare at the top
//	//  number of method stencil points ("Differenzsternpunkte")
//	int 			noStencilPoints = 3;
//	// dimension
	int 			dim = 3;
//	// increment
//	double 			d = 1e-6;
//	// turbulent perturbations
//    vector<double> 	utrp(noStencilPoints + 1);
//    vector<double> 	vtrp(noStencilPoints + 1);
//    vector<double> 	wtrp(noStencilPoints + 1);
//
//    // Loop over all stencil points
//	for (int i = 0; i <= noStencilPoints; i++)
//	{
		vector<double> 	p(dim);	// copy x to p, since x not modifiable
		p[0] = x[0];
		p[1] = x[1];
		p[2] = x[2];

//		#makeFieldInc
//		if ( i != 0 )
//		{
//		p[i-1] += d;
//		}

		// Set the turbulent velocities to zero for initialisation
		// TODO: #makeFieldInc was not included!!!
		utrp = 0;
		vtrp = 0;
		wtrp = 0;


	    // Loop over all wavenumbers
	    for (int m = 0; m <= nmodes-1; ++m)
		{
			kxi = kxio[m];
			kyi = kyio[m];
			kzi = kzio[m];

			sx = sxio[m];
			sy = syio[m];
			sz = szio[m];

			kx = kxi*wnr[m];
			ky = kyi*wnr[m];
			kz = kzi*wnr[m];
			rk = sqrt(pow(kx,2)+pow(ky,2)+pow(kz,2));

			// if the wavenumber, rk, is smaller than the largest wavenumber, then create fluctuations
			if (rk < wnrn)
			{
				arg = kx*p[0]+ky*p[1]+kz*p[2]+psi[m];
				tfunk = cos(arg);

				// modified von Karman spectrum
				e = amp/wnre*pow(wnr[m]/wnre,4)/pow(1+pow(wnr[m]/wnre,2),17./6.)*exp(-2*pow(wnr[m]/wnreta,2));

				utn = sqrt(e*pow(up,2)*dkn[m]);

				// turbulent perturbations
				utrp += 2*utn*tfunk*sx;
				vtrp += 2*utn*tfunk*sy;
				wtrp += 2*utn*tfunk*sz;

//				#makeFieldInc
//				// turbulent perturbations
//				utrp[i] += 2*utn*tfunk*sx;
//				vtrp[i] += 2*utn*tfunk*sy;
//				wtrp[i] += 2*utn*tfunk*sz;
			}  // end of if (int rk < wnrn)
		}  // end of for (int m = 0; m <= nmodes-1; ++m)
//	} #makeFieldInc

    //------------------------------------------------------------------------------------------------------------------------------------------------------
    // Blending function
    //------------------------------------------------------------------------------------------------------------------------------------------------------
    // ATTENTION: Currently works only for plane walls perpendicular to the y-axis.

//	// Linear blending for viscous sublayer
//	double uPlus_trp, vPlus_trp, wPlus_trp;
//	double yPlus = x[1] * ( ReTau / h );
//
//	if ( yPlus <= 5 )
//	{
//		uPlus_trp = yPlus;
//		utrp = uPlus_trp * ( ReTau * visc / h );
//		vtrp = uPlus_trp * ( ReTau * visc / h );
//		wtrp = uPlus_trp * ( ReTau * visc / h );
//	}

    // Distance to the nearest wall, TODO:not valid for the cell centres
    minDist = std::min(x[2], h-x[2]);

    // Manipulated hyperbolic function, freeStreamTurb prescribes a free-stream turbulence
    //	by preventing the blending function from dropping below the set value
    fBlend = std::max(0.5*(1 - tanh((minDist - delta)/blendDist)), freeStreamTurb);
    fBlend = 1;

//  #makeFieldInc
//	for (int i = 0; i <= noStencilPoints; i++)
//	{
//		utrp[i] *= fBlend;
//		vtrp[i] *= fBlend;
//		wtrp[i] *= fBlend;
//	}
//
//	//------------------------------------------------------------------------------------------------------------------------------------------------------
//	// Make perturbation field divergence free (incompressible)
//	//------------------------------------------------------------------------------------------------------------------------------------------------------
//	//TODO: declare at the top
//	vector<double>	gradU(dim), gradV(dim), gradW(dim);
//
//	// Calculate gradients
//	for (int i = 0; i < dim; i++)
//	{
//	// gradU
//		double f = utrp[0];
//		double f_d = utrp[i+1];
//		gradU[i] = ( (f_d - f) / d );
//	// gradV
//		f = vtrp[0];
//		f_d = vtrp[i+1];
//		gradV[i] = ( (f_d - f) / d );
//	// gradW
//		f = wtrp[0];
//		f_d = wtrp[i+1];
//		gradW[i] = ( (f_d - f) / d );
//	}
//
//	double utrp_inc, vtrp_inc, wtrp_inc; // incompressible turbulent perturbations
//
//	utrp_inc = gradW[1] - gradV[2]; // dW/dy - dV/dz
//	vtrp_inc = gradU[2] - gradW[0]; // dU/dz - dW/dx
//	wtrp_inc = gradV[0] - gradU[1]; // dV/dx - dU/dy

    double Uin = m_flow->getInletVelocity();
	double uPlus_in;
	double h_half = h/2;
	double yPlus; // ReTau = y*(ReTau/h_half);

	// yPlus for center line and upper & lower channel part
	if ( x[2] == h_half)
	{
		yPlus = ReTau;
	}
	else
	{
		yPlus = minDist * ( ReTau / h_half );
	}
	// mean inlet velocity profile towards walls
	if ( yPlus <= 5 )
	{
		uPlus_in = yPlus;
	}
	else if ( yPlus > 5 && yPlus < 30)
	{
		uPlus_in = -3.05 + 5. * log(yPlus);
	}
	else
	{
		uPlus_in = 1./0.4 * log(yPlus) + 5.2;
	}
	Uin = uPlus_in * ( ReTau * visc / h_half );

    // Vin & Win are assumed to be 0.
	if (component == 0)
	{
		return ( Uin + fBlend*utrp );
	}
	else if (component == 1)
	{
		return ( fBlend*vtrp );
		//return ( vtrp );
	}
	else // component == 2
	{
		return ( fBlend*wtrp );
		//return ( wtrp );
	}

    /*
    // TODO: add initialisation cases: InletBC_id
	if (component == 0) {
		// exponential law profile
		if (x(2) <= h/2)
			return ( m_flow->m_uBulk * std::pow(2*x(2)/h, 1./7.) );
		else
			return ( m_flow->m_uBulk * std::pow(2*(1-x(2)/h), 1./7.) );

		// parabolic profile
		return (- 4 * m_flow->m_uBulk * 1.5 *
				(x(2) - h) * x(2) / (h*h) );
		*/
} // end of InitialVelocity



//============================ MEAN ===========================

inline void TurbulentChannelFlow3D::mean(vector<vector<double> > &matrix, double &avg){

	int 				i = matrix.size(), j = matrix[0].size();
	vector<double> 		sum(i);

	for (int k = 0; k <= i-1; ++k){
		sum[k] = accumulate(matrix[k].begin(), matrix[k].end(), 0.0);
	}
	avg = accumulate(sum.begin(), sum.end(), 0.0) / (i*j);
}


//========================== ANGLES ===========================

inline void TurbulentChannelFlow3D::angles(int nmodes, int &iseed, int &iy, vector<int> &iv, vector<double> &fi,
		vector<double> &psi, vector<double> &alfa, vector<double> &teta){

vector<double> ang(nmodes);

randf_1(nmodes, 0, 2*M_PI, iseed, iy, iv, fi, iseed);
randf_1(nmodes, 0, 2*M_PI, iseed, iy, iv, psi, iseed);
randf_1(nmodes, 0, 2*M_PI, iseed, iy, iv, alfa, iseed);
randf_1(nmodes, 0, 1, iseed, iy, iv, ang, iseed);

for (int m = 0; m <= nmodes-1; ++m){
	teta[m] = acos(1-ang[m]/0.5);
}
}


//========================== RANDF_1 ===========================

inline void TurbulentChannelFlow3D::randf_1(int nmd, int alow, double ahigh, int idum, int &iy, vector<int> &iv,
		vector<double> &out, int &iseed){

double 	ran1 	= 0;

for (int j = 0; j <= nmd-1; ++j){
	randf_2(idum, iy, iv, ran1, iseed);
	idum = iseed;
	out[j] = (ran1*(ahigh-alow))+alow;
}
}


//========================== RANDF_2 ===========================

/* Generates a random number 0 < x < 1

   Minimal random number generator of Park and Miller with Bays-Durham shuffle and
   added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
	inline void timeNow(int &iseed);   the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
   alter idum between successive deviates in a sequence. RNMX should approximate the largest
   floating value that is less than 1.

   Park, S.K., and Miller, K.W. 1988, Communications of the ACM, vol. 31, pp. 1192-1201.
*/

inline void TurbulentChannelFlow3D::randf_2(int idum, int &iy, vector<int> &iv, double &ran1, int &iseed){

int			ia		= 16807;
double		im 		= 2147483647;
double		am		= 1/im;
int 		iq		= 127773;
int 		ir		= 2836;
int 		ntab	= 32;
double 		ndiv	= 1+(im-1)/ntab;
double 		eps		= 1.2e-7;
double		rnmx	= 1-eps;

int 			j, k;

//initial iseed (idum) is negative
if (idum <= 0 || iy == 0){
	idum = std::max(-idum, 1);
   	for (int j = ntab+8; j >=1; --j){
      		k = floor(idum/iq);
      		idum = ia*(idum-k*iq)-ir*k;
      		if (idum < 0){
         		idum = idum+im;
      		}
      		if (j <= ntab){
         		iv[j-1] = idum;
		}
	}
	iy = iv[0];
}

k 		= floor(idum/iq);
idum 	= ia*(idum-k*iq)-ir*k;

if  (idum <= 0){
	idum = idum+im;
}

j 		= floor(iy/ndiv);
iy 		= iv[j];
iv[j]	= idum;
iseed	= idum;
ran1	= std::min(am*iy,rnmx);
}

} /* namespace natrium */

