#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <iostream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor


//#define RE 1.0
#define H 0.3
#define UA 0.1
#define PI 3.141592653589793238463

// Sinus profile for top wall
template<typename T>
class SinusShapeDomain2D : public DomainFunctional2D {
public:
    SinusShapeDomain2D(T h_, T lx_, T b_): h(h_), lx(lx_), b(b_)
    { }
    virtual bool operator() (plint iX, plint iY) const {
        return iY > b * sin(2 * PI * iX / lx) + h;
    }
    virtual SinusShapeDomain2D<T>* clone() const {
        return new SinusShapeDomain2D<T>(*this);
    }
private:
    T h;
    T lx;
    T b;

};


void channelSetup (
        MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
        IncomprFlowParam<T> const& parameters,
        OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
        plint n, T b, T lx)
{
    // Create Velocity boundary conditions.
    //boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    // Specify the boundary velocity.
    Array<T,2> u (UA,0);
    Array<T,2> u0(0,0);


    const plint nx = lattice.getNx();

    lattice.periodicity().toggle(0, true); // Set periodic boundaries.

    // Moving bottom wall
    Box2D bottomWall(0, nx, 0, 0);
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottomWall );	// Linear condition on the bottom wall
    setBoundaryVelocity (lattice, bottomWall, u);

    // Sinus Profile
    defineDynamics(lattice, lattice.getBoundingBox(),
                   new SinusShapeDomain2D<T>(H * n, lx * n, b * n),
                   new BounceBack<T,DESCRIPTOR>);

    // Create the initial condition.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), 1.0, u0 );

    lattice.initialize();
}

void writeGifs(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
}

T computeH(plint iX, plint n, T b, T lx)
{
    return b * sin(2 * PI * iX / (lx * n)) + H;
}



T computeAverageULine(auto_ptr<MultiTensorField2D<T,2> > u, plint iX, plint n, T b, T lx)
{
    T result = 0;
    plint c = 0;
    for (plint iY = 0; iY < u->getNy(); iY++)
    {
        if(iY <= computeH(iX,n,b,lx) * n)
        {
            result += u->get(iX,iY)[0];
            c++;
        }
    }
    result = UA - result / c;
    return result;
}

T computeAverageU(auto_ptr<MultiTensorField2D<T,2> > u, plint n, T b, T lx)
{
    T result = 0;
    plint c = 0;
    for (plint iX = 0; iX < u->getNx(); iX++)
    {
        for (plint iY = 0; iY < u->getNy(); iY++)
        {
            if((iY <= computeH(iX,n,b,lx) * n))
            {
                result += u->get(iX,iY)[0];
                c++;
            }
        }
    }
    result = UA - result / c;
    return result;
}

T computeFlowFactor(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, T lx, T b, plint n, T sigma)
{
    T flowFactor = 0;

    for (plint iX = 0; iX < lattice.getNx(); iX++)
    {
        T averageU = computeAverageULine(computeVelocity(lattice),iX,n,b,lx);
        flowFactor += (2 * averageU / UA - 1.0) * computeH(iX,n,b,lx) / sigma;
    }
    return flowFactor / (lattice.getNx());
}

/**<  T computeFlowFactor(auto_ptr<MultiTensorField2D<T,2> > u)
{
    T sigma = B / sqrt(2.);
    T tmp = 0;

    for (plint iX = 0; iX < LX * N; iX++)
    {
        T averageU = computeAverageULine(u, iX);
        tmp += (2 * averageU / UA - 1.0) * computeH(iX) / sigma;
    }
    return tmp / LX;
} */

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    plint N;
	T RU, LX, B, RE;

	// Check the number of parameters
	if (argc < 6) {
		// Tell the user how to run the program
		pcerr << "Usage: " << argv[0] << " N Ru Lx b RE" << endl;
		/* "Usage messages" are a conventional way of telling the user
		* how to run a program if they enter the command incorrectly.
		*/
		return 1;
	}
	N = atoi(argv[1]);
	RU = atof(argv[2]);
	LX = atof(argv[3]);
    B = atof(argv[4]);
    RE = atof(argv[5]);

    // Use the class IncomprFlowParam to convert from
    //   dimensionless variables to lattice units, in the
    //   context of incompressible flows.
    IncomprFlowParam<T> parameters(
            RU,  // Reference velocity (the maximum velocity
                       //   in the Poiseuille profile) in lattice units.
            (T) RE,  // Reynolds number
            N,       // Resolution of the reference length (channel height).
            LX,        // Channel length in dimensionless variables
            H + B         // Channel height in dimensionless variables
    );


    // Create lattice
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
             parameters.getNx(), parameters.getNy(),
             new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );
    lattice.toggleInternalStatistics(false);
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    channelSetup(lattice, parameters, *boundaryCondition, N, B, LX);

    // Abort criterion
    T sigma = B / sqrt(2.);
    T err = 1;
    //auto_ptr<MultiTensorField2D<T,2> > uLast = computeVelocity(lattice);
    T lastFlowFactor = 100;

    pcout << "Reference velocity: " << RU << endl
          << "Reynolds number:    " << RE << endl
          << "Resolution:         " << N << endl
          << "Channel length:     " << LX << endl
          << "Channel height:     " << H+B << endl
          << "Sin-Amplitude:      " << B << endl;

    // Main loop over time iterations.
   while (err > 1e-6) {

        // Execute lattice Boltzmann iteration.
        for(plint i = 0; i < 100; i++)
            lattice.collideAndStream();

	lattice.toggleInternalStatistics(true);
        //T flowFactor = computeFlowFactor(lattice, LX, B, N, sigma);
        T u_average = computeAverage( *computeVelocityComponent(lattice,0));
        // correction due to averaging over whole rectangle
        u_average *= ceil((H+B)*N) -2;
        u_average /= ceil(H*N) -2;
        u_average = UA - u_average;

        T flowFactor = 2*u_average *H / (sigma * UA) - H / sigma;
        //pcout << flowFactor << " " << flowFactor2 << endl;

        err = abs(flowFactor - lastFlowFactor);
        lastFlowFactor = flowFactor;
        /*for (plint iX = 0; iX < u->getNx(); iX++)
        {
            for (plint iY = 0; iY < u->getNy(); iY++)
            {
                T tempUx = u->get(iX,iY)[0] - uLast->get(iX,iY)[0];
                T tempUy = u->get(iX,iY)[1] - uLast->get(iX,iY)[1];
                T tempErr = sqrt(tempUx*tempUx + tempUy*tempUy);
                if(tempErr > err) err = tempErr;
            }
        }*/


        //writeGifs(lattice, iT);
        //uLast = computeVelocity(lattice);

	lattice.toggleInternalStatistics(false);

    }

    //Flowfactor
    //T averageU = computeAverageU(computeVelocity(lattice));

    T flowFactor = computeFlowFactor(lattice, LX, B, N, sigma);
    T epsilon = (H - B) / (B / sqrt(2.));
    pcout << "Average u:          " << computeAverageU(computeVelocity(lattice),N,B,LX) << endl
          << "Flowfactor:         "  << flowFactor << endl
          << "Epsilon:            " << epsilon << endl
          << "Tau:                " << parameters.getTau() << endl;
    delete boundaryCondition;
}
