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
#include "deal.II/grid/grid_tools.h"
#include "deal.II/base/tensor.h"

#include "natrium/problemdescription/LinearBoundaryRhoU.h"
#include "natrium/problemdescription/PeriodicBoundary.h"
#include "natrium/problemdescription/NonlinearBoundaryZouHeRho.h"
#include "natrium/problemdescription/ConstantExternalForce.h"
#include "natrium/utilities/Math.h"

namespace natrium {

TurbulentChannelFlow3D::TurbulentChannelFlow3D(double viscosity, size_t refinementLevel,
		double ReTau, double ReCl, double u_bulk, double height, double length, double width,
		double orderOfFiniteElement, bool is_periodic) :
		ProblemDescription<3>(makeGrid(height, length, width), viscosity, height),
		m_refinementLevel(refinementLevel), m_ReTau(ReTau), m_ReCl(ReCl), m_uBulk(u_bulk),
		m_ofe(orderOfFiniteElement), m_height(height), m_length(length), m_width(width),
		m_maxUtrp(0.0), m_maxIncUtrp(0.0) {

	// Recommendations for CPU use
	cout << "-------------------------------------------------------------" << endl;
	cout << "**** Recommendations for CPU use ****" << endl;
	double noCube = ( width/height ) * ( length/height );
	double noGridPoints = pow( (orderOfFiniteElement + 1), 3 ) * pow(8, refinementLevel) * noCube;
	cout << "... Computation node details: " << endl;
	cout << "    - #CPU per node: 12 " << endl;
	cout << "    - memory per node: 4000 MB " << endl;
	cout << "... Predicted number of total grid points: " << noGridPoints << endl;
	cout << "... Recommended number of total grid points per node: 10e+6" << endl;
	cout << "... Recommended number of nodes: " << ceil(noGridPoints/10e+6) << endl;
	cout << "------------------------------------------------------------" << endl;

	/// apply boundary values
	setBoundaries(makeBoundaries(is_periodic));
	// apply initial values / analytic solution
	setInitialU(boost::make_shared<IncompressibleU>(this));

	if (is_periodic) {
		// add external force
		// turbulent flow
		double rho = 1;
		double h_half = 0.5 * height;
		double Fx = pow( ReTau * viscosity / h_half, 2) * rho / h_half;
		// laminar flow
		//double Fx = 8 * 1.5 * m_uBulk * viscosity / (height * height);
		cout << " >>>> Body force F = " << Fx << endl;
		dealii::Tensor<1, 3> F;
		F[0] = Fx;
		setExternalForce(
				boost::make_shared<ConstantExternalForce<3> >(F,
						SHIFTING_VELOCITY));
	}

	// refine global
	getMesh()->refine_global(refinementLevel);
	// transform grid to unstructured grid
	dealii::GridTools::transform(UnstructuredGridFunc(height), *getMesh());
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
	repetitions.at(0) = (int) length / height;
	repetitions.at(1) = (int) height;
	repetitions.at(2) = (int) width / height;
	dealii::GridGenerator::subdivided_hyper_rectangle(*mesh, repetitions,
			dealii::Point<3>(0.0, 0.0, 0.0),
			dealii::Point<3>(length, height, width), colorize);

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

	//cout << "**** Set boundary conditions ****" << endl;
	if (is_periodic) {
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh())); // TODO: damn, comment it!
		//cout << " > periodic: inlet/outlet" << endl;
		boundaries->addBoundary(
				boost::make_shared<PeriodicBoundary<3> >(4, 5, 2, getMesh()));
		//cout << " > periodic: back/front" << endl;
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(2, zeroVector));
		//cout << " > no-slip: top" << endl;
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(3, zeroVector));
		//cout << " > no-slip: bottom" << endl;

	} else {
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(0,
						boost::make_shared<dealii::ConstantFunction<3> >(1.0),
						boost::make_shared<TurbulentChannelFlow3D::InitialVelocity>(this)));
		boundaries->addBoundary(
				boost::make_shared<NonlinearBoundaryZouHeRho<3> >(1,
						boost::make_shared<dealii::ConstantFunction<3> >(-1.0), 1));
		boundaries->addBoundary(
						boost::make_shared<PeriodicBoundary<3> >(4, 5, 2, getMesh()));
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(2, zeroVector));
		boundaries->addBoundary(
				boost::make_shared<LinearBoundaryRhoU<3> >(3, zeroVector));
	}

	return boundaries;
}


double TurbulentChannelFlow3D::MeanVelocityProfile::value(const dealii::Point<3>& x,
		const unsigned int component) const{

	//return m_initialIncompressibleU.value(x, component);

	double visc 	= m_flow->getViscosity();
    double ReTau 	= m_flow->getFrictionReNumber();
	double ReCl 	= m_flow->getCenterLineReNumber();
	double height 	= m_flow->getCharacteristicLength();

    double Uin = ReCl * 2 * visc / height;
	double uPlus_in;

    double minDist = std::min(x[1], height - x[1]);
	double h_half = 0.5 * height; // half channel height

	// yPlus towards upper & lower channel wall
	double yPlus = ReTau * ( minDist / h_half );

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

    // mean velocities <Vin> & <Win> are assumed to be 0.
	if (component == 0)
	{
		return Uin;
		//return ( Uin + m_initialIncompressibleU.value(x, component) );
	}
	else // component 1 & 2
	{
		return 0;
		//return ( m_initialIncompressibleU.value(x, component) );
	}
}


double TurbulentChannelFlow3D::IncompressibleU::value(const dealii::Point<3>& x,
		const unsigned int component) const{

	//return m_initialU.value(x, component);

	//---------------------------------------------------------
	// Make perturbation field divergence free (incompressible)
	//---------------------------------------------------------
	const int		dim = 3; // TODO: applicable both to 3D and 2D
	vector<double> 	gradU(dim, 0.0), gradV(dim, 0.0), gradW(dim, 0.0);

	dealii::Point<dim>	x_h = x;

	// Calculate gradients
	for (int componentID = 0; componentID < dim; componentID++)
	{
		x_h[componentID] = x[componentID] + m_increment;
	// gradU
		double f = m_initialU.value(x, 0);
		double f_h = m_initialU.value(x_h, 0);
		gradU[componentID] = ( (f_h - f) / m_increment );
	// gradV
		f = m_initialU.value(x, 1);
		f_h = m_initialU.value(x_h, 1);
		gradV[componentID] = ( (f_h - f) / m_increment );
	// gradW
		f = m_initialU.value(x, 2); // z is the height coordinate!
		f_h = m_initialU.value(x_h, 2);
		gradW[componentID] = ( (f_h- f) / m_increment );

		x_h[componentID] = x[componentID];
	}

	// incompressible turbulent perturbations
	double utrp_inc = gradW[1] - gradV[2]; // dW/dy - dV/dz
	double vtrp_inc = gradU[2] - gradW[0]; // dU/dz - dW/dx
	double wtrp_inc = gradV[0] - gradU[1]; // dV/dx - dU/dy

	if ( utrp_inc > m_flow->m_maxIncUtrp )
	{
		m_flow->m_maxIncUtrp = utrp_inc;
		//cout << "maxIncUtrp = "<< m_flow->m_maxIncUtrp << endl;
	}

	if (component == 0)
	{
		return ( utrp_inc );
	}
	else if (component == 1)
	{
		return ( vtrp_inc );
	}
	else // component == 2
	{
		return ( wtrp_inc );
	}
}

double TurbulentChannelFlow3D::InitialVelocity::value(const dealii::Point<3>& x,
		const unsigned int component) const {

	//pout << "here" << endl;

	// Synthetic turbulence generation due to Davidson et al.
	// [1] Using isotropic...
	// [2]

	// turbulent kinetic energy to be set properly

	assert(component < 3);

	//------------------------------------------------------------------------------------------------------------------------------------------------------
	// Selectable variables
	// TODO: These variables have to be made selectable from a property file or similar

	double 	height = m_flow->getCharacteristicLength();				// characteristic flow length scale, here: full channel height
	double  visc = m_flow->getViscosity();						// kinematic viscosity
	//double	ReTau = m_flow->getFrictionReNumber();
	double	qm 					= 4.0;							// turbulent kinetic energy
	double	delta 				= 5./32 * height;	        			// inlet boundary layer thickness. Only if the boundary layer
																// at inlet is not fully developed
	//double 	uTau 				= 1/25.0;						// inlet shear velocity, suTau = Urms [Ref2]

	// blending function parameters
	double	blendDist 			= 0.1 * height;						// distance over which fBlend goes from 0 to 1
	double	freeStreamTurb 		= 0.1;							// parameter does not let fBlend drop below the prescribed value

	int 	nmodes 				= 32;							// number of Fourier modes
	double	wew1fct				= 2;							// ratio of ke and kmin (in wavenumber)

	// Constants
	double 	amp 				= 1.452762113;					// alpha in [Ref2]

	// Calculated
	double	sli 				= 0.1*delta;    				// length scale
	double	urms	 			= sqrt(2*qm/3);					// turbulent velocity scale
	double	epsm 				= pow(qm,1.5)/sli;				// dissipation rate

	double 	ofe					= m_flow->getOrderOfFiniteElement();
	double	inletCellLength		= height/( (ofe + 1) * pow(2, m_flow->getRefinementLevel()) );
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

	// Set the turbulent velocities to zero for initialisation
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
		rk = sqrt( pow(kx, 2) + pow(ky, 2) + pow(kz, 2) );

		// if the wavenumber, rk, is smaller than the largest wavenumber, then create fluctuations
		if (rk < wnrn)
		{
			arg = kx*x[0] + ky*x[1] + kz*x[2] + psi[m];
			tfunk = cos(arg);

			// modified von Karman spectrum
			e = amp/wnre * pow( wnr[m]/wnre, 4 ) / pow( 1 + pow( wnr[m]/wnre, 2 ), 17./6. )
					* exp( -2 * pow( wnr[m]/wnreta , 2 ) );

			utn = urms * sqrt( e * dkn[m] );

			// turbulent perturbations
			utrp += 2 * utn * tfunk * sx;
			vtrp += 2 * utn * tfunk * sy;
			wtrp += 2 * utn * tfunk * sz;
		}  // end of if (int rk < wnrn)
	}  // end of for (int m = 0; m <= nmodes-1; ++m)

	// seek max(utrp);
	if ( utrp > m_flow->m_maxUtrp )
	{
		m_flow->m_maxUtrp = utrp;
		//cout << "maxUtrp = "<< m_flow->m_maxUtrp << endl;
	}


    //------------------------------------------------------------------------------------------------------------------------------------------------------
    // Blending function
    //------------------------------------------------------------------------------------------------------------------------------------------------------
    // ATTENTION: Currently works only for plane walls perpendicular to the y-axis.

    // Distance to the nearest wall, TODO:not valid for the cell centres
    minDist = std::min(x[1], height - x[1]);

/*
    // Linear blending for viscous sublayer
	double h_half = 0.5 * height; // half channel height
	double yPlus = ReTau * ( minDist / h_half );

	// TODO: min(1,1) = ?;

	if ( yPlus <= 5 )
	{
		double uPlus_trp = yPlus;
		utrp = uPlus_trp * ( ReTau * visc / h_half );
		vtrp = uPlus_trp * ( ReTau * visc / h_half );
		wtrp = uPlus_trp * ( ReTau * visc / h_half );
	}
*/

    // Manipulated hyperbolic function, freeStreamTurb prescribes a free-stream turbulence
    //	by preventing the blending function from dropping below the set value
    fBlend = std::max( 0.5 * ( 1 - tanh( (minDist - delta)/blendDist ) ), freeStreamTurb);
    //fBlend = 1;

    // DEBUG: for mean velocity profile check
    //return 0;

    // Vin & Win are assumed to be 0.
	if (component == 0)
	{
		return ( fBlend*utrp );
	}
	else if (component == 1)
	{
		return ( fBlend*vtrp );
	}
	else // component == 2
	{
		return ( fBlend*wtrp );
	}

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

