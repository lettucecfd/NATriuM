/*
 * MixingLayer3D.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */

#include "MixingLayer3D.h"

#include "deal.II/grid/grid_generator.h"
//#include "deal.II/grid/tria_accessor.h"
//#include "deal.II/grid/tria_iterator.h"
//#include "deal.II/grid/grid_out.h"
//#include "deal.II/fe/fe_series.h"
//#include "../postprocessing/FFT.cpp"
//#include "../postprocessing/fftw3.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"
#include <random>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>
//#include <complex>
//#include <vector>
//#include <iomanip>

using namespace std;

double shearlayerthickness = 0.093;
int kmax = 48; // [1] C. Pantano and S. Sarkar, “A study of compressibility effects in the high-speed turbulent shear layer using direct simulation,” J. Fluid Mech., vol. 451, pp. 329–371, Jan. 2002, doi: 10.1017/S0022112001006978.


namespace natrium {

MixingLayer3D::MixingLayer3D(double viscosity,
                             size_t refinementLevel, double U) :
        ProblemDescription<3>(makeGrid(), viscosity, 1), m_U(U), m_refinementLevel(refinementLevel) {
    /// apply boundary values
    setBoundaries(makeBoundaries());
    // apply analytic solution
    this->setInitialU(boost::make_shared<InitialVelocity>(this));
    this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

MixingLayer3D::~MixingLayer3D() = default;

double MixingLayer3D::InitialVelocity::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component < 3);
    double rand_u = InterpolateVelocities(x(0), x(1), x(2), component);
    rand_u *= exp(-pow((x(1))/(2*shearlayerthickness),2));
    if (component == 0) {
        return tanh(-x(1)/(2*shearlayerthickness)) + rand_u;
    } else {
        return rand_u;
    }
}

MixingLayer3D::InitialVelocity::InitialVelocity(natrium::MixingLayer3D *flow) : m_flow(flow) {
    vector<double> x, y, z;
    double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz;
    float lx, ly, lz;
    int npoints = 20;
    lx = 1720;
    ly = 387;
    lz = 172;
    xmax = lx * shearlayerthickness / 2;
    xmin = -xmax;
    ymax = ly * shearlayerthickness / 2;
    ymin = -ymax;
    zmax = lz * shearlayerthickness / 2;
    zmin = -zmax;
    dx = min({(xmax-xmin)/npoints, (ymax-ymin)/npoints, (zmax-zmin)/npoints});
    dy = dx;
    dz = dx;
    nx = ceil((xmax-xmin)/dx);
    ny = ceil((ymax-ymin)/dy);
    nz = ceil((zmax-zmin)/dz);

    // create linspaces x, y, z
    double linvalue;
    linvalue = xmin; for (int i = 0; i < nx; i++) { x.push_back(linvalue); linvalue += dx; }
    linvalue = ymin; for (int i = 0; i < ny; i++) { y.push_back(linvalue); linvalue += dy; }
    linvalue = zmin; for (int i = 0; i < nz; i++) { z.push_back(linvalue); linvalue += dz; }
    xvec.insert(xvec.end(), x.begin(), x.end());
    yvec.insert(yvec.end(), y.begin(), y.end());
    zvec.insert(zvec.end(), z.begin(), z.end());

    auto xbounds = std::minmax_element(x.begin(), x.end());
    minx = *xbounds.first; maxx = *xbounds.second;
    auto ybounds = std::minmax_element(y.begin(), y.end());
    miny = *ybounds.first; maxy = *ybounds.second;
    auto zbounds = std::minmax_element(z.begin(), z.end());
    minz = *zbounds.first; maxz = *zbounds.second;

    // warm up randomness
    static int sseq[ std::mt19937::state_size ] ;
    const static bool once = ( std::srand( std::time(nullptr)), // for true randomness: std::time(nullptr)
            std::generate( std::begin(sseq), std::end(sseq), std::rand ),
            true ) ;
    static std::seed_seq seed_seq( std::begin(sseq), std::end(sseq) ) ;
    static std::mt19937 twister(seed_seq) ;
    // random generator in [-1,1]
    static uniform_real_distribution<double> distr(-1.0, 1.0) ; // +- Velocity

    double k0 = 23.66 * shearlayerthickness; // peak wave number

    // Fill randomPsi with random values
//    vector< vector< vector< vector<double> > > > randomPsi;
    randomPsi.reserve(3);
//    randomPsi.reserve(3);
    for (int dir = 0; dir < 3; dir++) { vector< vector< vector<double> > > tmpdir;
        for (int xi = 0; xi < nx; xi++) {
            vector<vector<double> > tmpi;
            for (int yi = 0; yi < ny; yi++) {
                vector<double> tmpj;
                for (int zi = 0; zi < nz; zi++) {
                    double tmpk;
                    tmpk = distr(twister); // * shearlayerthickness;
                    tmpj.push_back(tmpk);
                } tmpi.push_back(tmpj);
            } tmpdir.push_back(tmpi);
        }
        // perform dft on randomPsi
        vector<vector<vector<complex<double>>>> psi_hat(nx, vector<vector<complex<double>>>(ny, vector<complex<double>>(nz))); // coordinates are on last three dimensions
        // dft along x
        for (int xi = 0; xi < nx; xi++) { vector<double> tmpfield;
            for (int yi = 0; yi < ny; yi++) {
                for (int zi = 0; zi < nz; zi++) { double tmp;
                    tmp = tmpdir[xi][yi][zi];
                    tmpfield.push_back(tmp);
                } }
            vector<complex<double>> dft_x = Fourier1D(tmpfield);
            for (int yi = 0; yi < ny; yi++) {
                for (int zi = 0; zi < nz; zi++) { complex<double> tmp;
                    tmp = dft_x[yi * nz + zi];
                    psi_hat[xi][yi][zi] += tmp;
                } }
        }
        // dft along y
        for (int yi = 0; yi < ny; yi++) { vector<double> tmpfield;
            for (int xi = 0; xi < nx; xi++) {
                for (int zi = 0; zi < nz; zi++) { double tmp;
                    tmp = tmpdir[xi][yi][zi];
                    tmpfield.push_back(tmp);
                } }
            vector<complex<double>> dft_y = Fourier1D(tmpfield);
            for (int xi = 0; xi < nx; xi++) {
                for (int zi = 0; zi < nz; zi++) { complex<double> tmp;
                    tmp = dft_y[xi * nz + zi];
                    psi_hat[xi][yi][zi] += tmp;
                } }
        }
        // dft along z
        for (int zi = 0; zi < nz; zi++) { vector<double> tmpfield;
            for (int yi = 0; yi < ny; yi++) {
                for (int xi = 0; xi < nx; xi++) { double tmp;
                    tmp = tmpdir[xi][yi][zi];
                    tmpfield.push_back(tmp);
                } }
            vector<complex<double>> dft_z = Fourier1D(tmpfield);
            for (int yi = 0; yi < ny; yi++) {
                for (int xi = 0; xi < nx; xi++) { complex<double> tmp;
                    tmp = dft_z[xi * ny + yi];
                    psi_hat[xi][yi][zi] += tmp;
                } }
        }
        // multiply in fourier space
        for (int xi = 0; xi < nx; xi++) {
            for (int yi = 0; yi < ny; yi++) {
                for (int zi = 0; zi < nz; zi++) { double k_abs;
                    k_abs = sqrt(xi*xi + yi*yi + zi*zi);
                    psi_hat[xi][yi][zi] *= exp(-2*k_abs/k0);
                }
            }
        }
        // perform inverse dft on psi_hat (directionally)
        // idft along x
        for (int xi = 0; xi < nx; xi++) {
            for (int yi = 0; yi < ny; yi++) {
                tmpdir[xi][yi] = InverseFourier1D(psi_hat[xi][yi]); // = idft_x
            }
        }
        // idft along y and add to idft_x
        for (int yi = 0; yi < ny; yi++) {
            for (int xi = 0; xi < nx; xi++) { vector<complex<double>> tmpfield(nz);
                for (int zi = 0; zi < nz; zi++) {
                    tmpfield[zi] = psi_hat[xi][yi][zi];
                }
                vector<double> idft_y = InverseFourier1D(tmpfield);
                for (int zi = 0; zi < nz; zi++) {
                    tmpdir[xi][yi][zi] += idft_y[zi];
                }
            }
        }
        // idft along z and add to idft_x and idft_y
        for (int zi = 0; zi < nz; zi++) {
            for (int xi = 0; xi < nx; xi++) { vector<complex<double>> tmpfield(ny);
                for (int yi = 0; yi < ny; yi++) {
                    tmpfield[yi] = psi_hat[xi][yi][zi];
                }
                vector<double> idft_z = InverseFourier1D(tmpfield);
                for (int yi = 0; yi < ny; yi++) {
                    tmpdir[xi][yi][zi] += idft_z[yi];
                }
            }
        }
        randomPsi.push_back(tmpdir);
    } // so: randomPsi = {psix, psiy, psiz} ;

    // calculate gradient using central difference scheme
    vector< vector< vector< vector< vector<double> > > > > gradient(3, vector<vector<vector<vector<double>>>>(3, vector<vector<vector<double>>>(nx, vector<vector<double>>(ny, vector<double>(nz))))); // coordinates are on last three dimensions
    int il, jl, kl, iu, ju, ku;
    for (int dir_psi = 0; dir_psi < 3; dir_psi++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    if (i==0) {il = 1;} else il = i;
                    if (j==0) {jl = 1;} else jl = j;
                    if (k==0) {kl = 1;} else kl = k;
                    if (i==nx-1) {iu = nx-2;} else iu = i;
                    if (j==ny-1) {ju = ny-2;} else ju = j;
                    if (k==nz-1) {ku = nz-2;} else ku = k;
                    gradient[dir_psi][0][i][j][k] = (randomPsi[dir_psi][iu+1][j][k] - randomPsi[dir_psi][il-1][j][k]) / dx;
                    gradient[dir_psi][1][i][j][k] = (randomPsi[dir_psi][i][ju+1][k] - randomPsi[dir_psi][i][jl-1][k]) / dy;
                    gradient[dir_psi][2][i][j][k] = (randomPsi[dir_psi][i][j][ku+1] - randomPsi[dir_psi][i][j][kl-1]) / dz;
                }}}} // so: gradient = {gradient_psix, gradient_psiy, gradient_psiz} and gradient_psin = { dpsin/dx, dpsin/dy, dpsin/dz }

    // calculate curl using gradient values
    int m, n; // for indices of cross-product
    double tmp;
    for (int dir_curl = 0; dir_curl < 3; dir_curl++) {
        vector<vector<vector<double> > > tmpdir;
        if (dir_curl == 0) {
            m = 2;
            n = 1;
        }
        else if (dir_curl == 2) {
            m = 1;
            n = 0;
        }
        else {
            m = dir_curl - 1;
            n = dir_curl + 1;
        }
        for (int i = 0; i < nx; i++) {
            vector<vector<double> > tmpi;
            for (int j = 0; j < ny; j++) {
                vector<double> tmpj;
                for (int k = 0; k < nz; k++) {
                    tmp = (gradient[m][n][i][j][k] - gradient[n][m][i][j][k]) / 4;
                    tmpj.push_back(tmp);
                }
                tmpi.push_back(tmpj);
            }
            tmpdir.push_back(tmpi);
        }
        curlOfPsi.push_back(tmpdir);
    } // so: curlOfPsi = {ux, uy, uz}
}

vector<complex<double>> MixingLayer3D::InitialVelocity::Fourier1D(const vector<double> &in) {
    int n = in.size();
    vector<complex<double>> out(n);
    for (int k = 0; k < kmax; ++k) {
        out[k] = {0.0, 0.0}; // Initialize each component to zero
        for (int t = 0; t < n; ++t) {
            double angle = 2.0 * M_PI * k * t / n;
            complex<double> exp_term(std::cos(angle), std::sin(angle));
            out[k] += in[t] * exp_term; //std::real(in[k] * exp_term);
        }
        out[k] /= n;
    }
    return out;
}

vector<double> MixingLayer3D::InitialVelocity::InverseFourier1D(const vector<complex<double>> &in) {
    int n = in.size();
    vector<double> out(n);
    for (int t = 0; t < n; ++t) {
        out[t] = 0.0; // Initialize each component to zero
        for (int k = 0; k < kmax; ++k) {
            double angle = 2.0 * M_PI * k * t / n;
            complex<double> exp_term(std::cos(angle), std::sin(angle));
            out[t] += std::real(in[k] * exp_term);
        }
        out[t] /= n; // Scale the result
    }
    return out;
}

double MixingLayer3D::InitialVelocity::InterpolateVelocities(double xq, double yq, double zq, const unsigned int dim) const
{
    /*
     * Assumes that all abscissa are monotonically increasing values
     */
    xq = std::max(minx, std::min(xq, maxx));
    yq = std::max(miny, std::min(yq, maxy));
    zq = std::max(minz, std::min(zq, maxz));

    auto xupper = std::upper_bound(xvec.cbegin(), xvec.cend(), xq);
    int x1 = (xupper == xvec.cend()) ? xupper - xvec.cbegin() - 1 : xupper - xvec.cbegin();
    int x0 = x1 - 1;

    auto yupper = std::upper_bound(yvec.cbegin(), yvec.cend(), yq);
    int y1 = (yupper == yvec.cend()) ? yupper - yvec.cbegin() - 1 : yupper - yvec.cbegin();
    auto y0 = y1 - 1;

    auto zupper = std::upper_bound(zvec.cbegin(), zvec.cend(), zq);
    int z1 = (zupper == zvec.cend()) ? zupper - zvec.cbegin() - 1 : zupper - zvec.cbegin();
    auto z0 = z1 - 1;

    double xd = (xq - xvec[x0])/(xvec[x1] - xvec[x0]);
    double yd = (yq - yvec[y0])/(yvec[y1] - yvec[y0]);
    double zd = (zq - zvec[z0])/(zvec[z1] - zvec[z0]);

    double c000 = curlOfPsi[dim][x0][y0][z0];
    double c010 = curlOfPsi[dim][x0][y1][z0];
    double c100 = curlOfPsi[dim][x1][y0][z0];
    double c110 = curlOfPsi[dim][x1][y1][z0];

    double c001 = curlOfPsi[dim][x0][y0][z1];
    double c011 = curlOfPsi[dim][x0][y1][z1];
    double c101 = curlOfPsi[dim][x1][y0][z1];
    double c111 = curlOfPsi[dim][x1][y1][z1];

    double c00 = c000*(1 - xd) + c100*xd;
    double c01 = c001*(1 - xd) + c101*xd;
    double c10 = c010*(1 - xd) + c110*xd;
    double c11 = c011*(1 - xd) + c111*xd;

    double c0 = c00*(1 - yd) + c10*yd;
    double c1 = c01*(1 - yd) + c11*yd;

    double c = c0*(1 - zd) + c1*zd;
    return c;
}

double MixingLayer3D::InitialDensity::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component == 0);
    return 1.0;// + p / (m_flow->m_cs * m_flow->m_cs);
}
double MixingLayer3D::InitialTemperature::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component == 0);
    return 1.0;
}

/**
 * @short create triangulation for Compressible Mixing Layer flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > MixingLayer3D::makeGrid() {
    //Creation of the principal domain
    boost::shared_ptr<Mesh<3> > cube = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
    double lx = 1720 * shearlayerthickness / 2;
    double ly = 387 * shearlayerthickness / 2;
    double lz = 172 * shearlayerthickness / 2;
    dealii::Point<3> corner1(-lx, -ly, -lz);
    dealii::Point<3> corner2(lx, ly, lz);
    vector<unsigned int> rep;
    rep.push_back(1);
    rep.push_back(1);
    rep.push_back(1);
    dealii::GridGenerator::subdivided_hyper_rectangle(*cube, rep, corner1, corner2, true);
    return cube;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3> > MixingLayer3D::makeBoundaries() {
    // make boundary description
    boost::shared_ptr<BoundaryCollection<3> > boundaries = boost::make_shared<
            BoundaryCollection<3> >();

    // velocity vector moving forward
    dealii::Vector<double> plusVector(3);
    plusVector[0]=1.0;
    plusVector[1]=0.0;
    plusVector[2]=0.0;

    // velocity vector moving backward
    dealii::Vector<double> minusVector(3);
    minusVector[0]=-1.0;
    minusVector[1]=0.0;
    minusVector[2]=0.0;

    // set boundaries on top and bottom to move forward / backward
    boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3> >(2, plusVector));
    boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3> >(3, minusVector));

    // set a boundary between 0 and 1, and 4 and 5, with direction 0 (x) and 2 (z), respectively
    boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(0, 1, 0, getMesh()));
    boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3> >(4, 5, 2, getMesh()));

    // Get the triangulation object (which belongs to the parent class).
    boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();
    return boundaries;
}

} /* namespace natrium */