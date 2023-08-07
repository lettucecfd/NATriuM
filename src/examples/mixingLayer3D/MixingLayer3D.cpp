/*
 * MixingLayer3D.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */

#include "MixingLayer3D.h"

#include "deal.II/grid/grid_generator.h"
#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"
#include <random>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>
using namespace std;

double shearlayerthickness = 0.093;
double k0 = 23.66 * shearlayerthickness; // peak wave number
//int n = 5;
//int kmax = pow(2, n); // [1] C. Pantano and S. Sarkar, “A study of compressibility effects in the high-speed turbulent shear layer using direct simulation,” J. Fluid Mech., vol. 451, pp. 329–371, Jan. 2002, doi: 10.1017/S0022112001006978.
// kmax = 32
//int npoints = 32; // number of points in shortest axis of velocity field (lz, presumably)
//bool print = false;
//bool recalculate_psi = true;

namespace natrium {

MixingLayer3D::MixingLayer3D(double viscosity, size_t refinementLevel, bool squash, bool print, bool recalculate, double U) :
    ProblemDescription<3>(makeGrid(), viscosity, 1), m_squash(squash),
    m_U(U), m_refinementLevel(refinementLevel) {
    if (m_refinementLevel > 4) { m_print = false; }
    /// apply boundary values
    setBoundaries(makeBoundaries());
    // apply analytic solution
    this->setInitialU(boost::make_shared<InitialVelocity>(this, print, recalculate));
    this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

MixingLayer3D::~MixingLayer3D() = default;

double MixingLayer3D::InitialVelocity::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component < 3);
    double rand_u = InterpolateVelocities(x(0), x(1), x(2), component);
    rand_u *= exp(-pow((x(1))/(2 * shearlayerthickness), 2));
    if (component == 0) {
        return tanh(-x(1)/(2 * shearlayerthickness)) + rand_u;
    } else {
        return rand_u;
    }
}

MixingLayer3D::InitialVelocity::InitialVelocity(natrium::MixingLayer3D *flow, bool print, bool recalculate) :
m_flow(flow), m_print(print), m_recalculate(recalculate) {
    int kmax = 48;//pow(2, flow->m_refinementLevel);
    kxmax = kmax; // maybe lower kxmax and k3 max?
    kymax = kmax;
    kzmax = kmax;
    vector<double> x, y, z;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    lx = 1720 * shearlayerthickness;
    ly = 387 * shearlayerthickness;
    lz = 172 * shearlayerthickness;
    xmax = lx / 2;// * 1.1;
    xmin = -xmax;
    ymax = ly / 2;// * 1.1;
    ymin = -ymax;
    zmax = lz / 2;// * 1.1;
    zmin = -zmax;
    nx = pow(2, flow->m_refinementLevel);
    ny = nx;
    nz = nx;
    kxmax = min(nx, kxmax);
    kymax = min(ny, kymax);
    kzmax = min(nz, kzmax);
    dx = lx / nx;
    dy = ly / ny;
    dz = lz / nz;

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

    if (m_recalculate) {
        cout << "Recalculating random velocity." << endl;
        //// warm up randomness
        static int sseq[ std::mt19937::state_size ] ;
        const static bool once = ( std::srand( std::time(nullptr)), // for true randomness: std::time(nullptr)
                std::generate( std::begin(sseq), std::end(sseq), std::rand ), true ) ;
        static std::seed_seq seed_seq( std::begin(sseq), std::end(sseq) ) ;
        static std::mt19937 twister(seed_seq) ;
        // random generator in [-1,1]
        static uniform_real_distribution<double> distr(-1., 1.) ; // +- Velocity

        randomPsi.reserve(3);
        for (int dir = 0; dir < 3; dir++) { vector< vector< vector<double> > > tmpdir;
            //// Fill randomPsi with random values
            cout << "Creating random velocity vector potential in direction " << dir << "." << endl;
            for (int xi = 0; xi < nx; xi++) { vector<vector<double> > tmpi;
                for (int yi = 0; yi < ny; yi++) { vector<double> tmpj;
                    for (int zi = 0; zi < nz; zi++) { double tmpk;
                        tmpk = distr(twister);
                        tmpj.push_back(tmpk);
                    } tmpi.push_back(tmpj);
                } tmpdir.push_back(tmpi);
            }
            if (m_print) {
                cout << "Printing random velocity vector potential in direction " << dir << "." << endl;
                ofstream rand_file("random_psi.txt");
                rand_file << "axis_" << to_string(dir) << endl;
                for (int i = 0; i < nx; i++) { rand_file << "[";
                    for (int j = 0; j < ny; j++) { rand_file << "[";
                        for (int k = 0; k < nz; k++) {
                            rand_file << tmpdir[i][j][k] << " ";
                        } rand_file << "]";
                    } rand_file << "],";
                } rand_file << endl;
            }
            cout << "Fourier transforming random velocity vector potential in direction " << dir << "." << endl;
            //// perform dft on randomPsi
            vector< vector< vector<complex<double>>>> psi_hat = Fourier3D(tmpdir);
            if (m_print) {
                cout << "Printing Fourier transform of random velocity vector potential in direction " << dir << "." << endl;
                ofstream dft_file("psi_hat.txt");
                dft_file << "axis_" << to_string(dir) << endl;
                for (int i = 0; i < nx; i++) { dft_file << "[";
                    for (int j = 0; j < ny; j++) { dft_file << "[";
                        for (int k = 0; k < nz; k++) {
                            dft_file << psi_hat[i][j][k] << " ";
                        } dft_file << "]";
                    } dft_file << "],";
                } dft_file << endl;
            }
            //// multiply in spectral space
            cout << "Scaling Fourier transformed velocity vector potential in direction " << dir << "." << endl;
            double k0 = sqrt(kxmax*kxmax+kymax*kymax+kzmax*kzmax);
            for (int kx = 0; kx < int(kxmax/2); kx++) {
                for (int ky = 0; ky < int(kymax/2); ky++) {
                    for (int kz = 0; kz < int(kzmax/2); kz++) {
                        double k_abs;
                        k_abs = sqrt(kx * kx + ky * ky + kz * kz);
                        psi_hat[kx][ky][kz] *= exp(-2 * k_abs / k0); // Sarkar:  (k_abs/k0)*(k_abs/k0)*(k_abs/k0)*(k_abs/k0)*exp(-2 * k_abs / k0); // holger:
                        psi_hat[kxmax-1 - kx][kymax-1 - ky][kzmax-1 - kz] *= exp(-2 * k_abs / k0);
                    } } }
            psi_hat[0][0][0] = 0;
            if (m_print) {
                cout << "Printing scaled Fourier transformed velocity vector potential in direction " << dir << "." << endl;
                ofstream scaling_file("psi_hat_scaled.txt");
                scaling_file << "axis_" << to_string(dir) << endl;
                for (int i = 0; i < nx; i++) { scaling_file << "[";
                    for (int j = 0; j < ny; j++) { scaling_file << "[";
                        for (int k = 0; k < nz; k++) {
                            scaling_file << psi_hat[i][j][k] << " ";
                        } scaling_file << "]";
                    } scaling_file << "],";
                } scaling_file << endl;
            }
            //// perform inverse dft on psi_hat (directionally)
            cout << "Inversely Fourier transforming velocity vector potential in direction " << dir << "." << endl;
            vector<vector<vector<double>>> psi_i = InverseFourier3D(psi_hat);
            if (m_print) {
                cout << "Printing inverse Fourier transform in direction " << dir << "." << endl;
                ofstream idft_file("psi_idft.txt");
                idft_file << "axis_" << to_string(dir) << endl;
                for (int i = 0; i < nx; i++) { idft_file << "[";
                    for (int j = 0; j < ny; j++) { idft_file << "[";
                        for (int k = 0; k < nz; k++) {
                            idft_file << psi_i[i][j][k] << " ";
                        } idft_file << "]";
                    } idft_file << "],";
                } idft_file << endl;
            }
            // add tmpdir (psix, psiy, or psiz) to randomPsi
            randomPsi.push_back(psi_i);
        } // so: randomPsi = {psix, psiy, psiz} ;

        // calculate gradient using central difference scheme
        cout << "Calculating gradients of velocity vector potential." << endl;
        vector<vector<vector<vector<vector<double>>>>> gradient(3,
            vector<vector<vector<vector<double>>>>(3,
                vector<vector<vector<double>>>(nx,
                    vector<vector<double>>(ny,
                        vector<double>(nz))))); // coordinates are on last three dimensions
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
        if (m_print) {
            cout << "Printing gradients of velocity vector potential." << endl;
            ofstream grad_file("psi_grad.txt");
            for (int dir_psi = 0; dir_psi < 3; dir_psi++) { grad_file << "direction: " << to_string(dir_psi) << endl << "[";
                for (int dir_grad = 0; dir_grad < 3; dir_grad++) { grad_file << "[";
                    for (int i = 0; i < nx; i++) { grad_file << "[";
                        for (int j = 0; j < ny; j++) { grad_file << "[";
                            for (int k = 0; k < nz; k++) {
                                grad_file << gradient[dir_psi][dir_grad][i][j][k] << " ";
                            } grad_file << "]";
                        } grad_file << "],";
                    } grad_file << "];" << endl;
                } grad_file << "];" << endl;
            }
        }
        // calculate curl using gradient values
        cout << "Calculating curl of velocity vector potential." << endl;
        int m, n; // for indices of cross-product
        for (int dir_curl = 0; dir_curl < 3; dir_curl++) {
            vector<vector<vector<double> > > tmpdir;
            if (dir_curl == 0) { m = 2; n = 1; }
            else if (dir_curl == 2) { m = 1; n = 0; }
            else { m = dir_curl - 1; n = dir_curl + 1; }
            for (int i = 0; i < nx; i++) { vector<vector<double> > tmpi;
                for (int j = 0; j < ny; j++) { vector<double> tmpj;
                    for (int k = 0; k < nz; k++) { double tmpk;
                        tmpk = (gradient[m][n][i][j][k] - gradient[n][m][i][j][k]) / 4;
                        tmpj.push_back(tmpk);
                    } tmpi.push_back(tmpj);
                } tmpdir.push_back(tmpi);
            } curlOfPsi.push_back(tmpdir);
        } // so: curlOfPsi = {ux, uy, uz}
        if (m_print) {
            cout << "Printing curl of velocity vector potential." << endl;
            ofstream u_file("random_u.txt");
            for (int dir_psi = 0; dir_psi < 3; dir_psi++) { //u_file << "[";
                for (int i = 0; i < nx; i++) { //u_file << "[";
                    for (int j = 0; j < ny; j++) { //u_file << "[";
                        for (int k = 0; k < nz; k++) {
                            u_file << ' ' << curlOfPsi[dir_psi][i][j][k];
                        } if (j<nz-1) { u_file << ','; }
                    } if (i<nz-1) { u_file << ';'; }
                } u_file << endl;
            }
        }
    }
    else { // TODO: Read data from random file
//        int dir_psi = 0, i = 0, j = 0, k = 0;
        ifstream file("random_u.txt");
        string line_dir;
        while (getline(file, line_dir)) {
            stringstream linestream(line_dir);
            string cell_dir;
            vector<vector<std::vector<double>>> tmpdir;
            while(std::getline(linestream, cell_dir, ';')) {
                string cell_i;
                vector<vector<double>> tmpi;
                while(getline(linestream,cell_i,',')) {
                    string cell_j;
                    vector<double> tmpj;
                    while(getline(linestream, cell_j, ' ')) {
                        string cell_k;
                        double tmpk = stod(cell_k);
                        tmpj.push_back(tmpk);
                    } tmpi.push_back(tmpj);
                } tmpdir.push_back(tmpi);
            } curlOfPsi.push_back(tmpdir);
        }
    }
}

vector<vector<vector<std::complex<double>>>>
MixingLayer3D::InitialVelocity::Fourier3D(const vector<vector<vector<double>>> &in) {
    vector<vector<vector<complex<double>>>> out(
            kxmax, vector<vector<complex<double>>>(
                    kymax, vector<complex<double>>(
                            kzmax, {0.0, 0.0}))); // coordinates are on last three dimensions
    //double xi, yi, zi;
    complex<double> omeg1, omeg2, omeg3;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
                for (int kx = 0; kx < kxmax; kx++) { omeg1 = {0.0, -2 * M_PI * ix/nx * kx};
                    for (int ky = 0; ky < kymax; ky++) { omeg2 = {0.0, -2 * M_PI * iy/ny * ky};
                        for (int kz = 0; kz < kzmax; kz++) { omeg3 = {0.0, -2 * M_PI * iz/nz * kz};
                            out[kx][ky][kz] += in[ix][iy][iz] * exp(omeg1 + omeg2 + omeg3);
    }}}}}}
    return out;
}

vector<vector<vector<double>>> MixingLayer3D::InitialVelocity::InverseFourier3D(
        const vector<vector<vector<std::complex<double>>>> &in) {
    vector<vector<vector<double>>> out(nx,
        vector<vector<double>>(ny,
            vector<double>(nz, 0.0)));
    //double xi, yi, zi;
    complex<double> omeg1, omeg2, omeg3, tmp;
    double N = (nx * ny * nz);
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
                for (int kx = 0; kx < kxmax; kx++) { omeg1 = {0.0, 2 * M_PI * ix/nx * kx};
                    for (int ky = 0; ky < kymax; ky++) { omeg2 = {0.0, 2 * M_PI * iy/ny * ky};
                        for (int kz = 0; kz < kzmax; kz++) { omeg3 = {0.0, 2 * M_PI * iz/nz * kz};
                            tmp = in[kx][ky][kz] * exp(omeg1 + omeg2 + omeg3);
                            out[ix][iy][iz] += real(tmp) / N;
    }}}}}}
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
    boost::shared_ptr<Mesh<3> > rect = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
    double lx = 1720 * shearlayerthickness / 2;
    double ly = 387 * shearlayerthickness / 2;
    double lz = 172 * shearlayerthickness / 2;
    dealii::Point<3> corner1(-lx, -ly, -lz);
    dealii::Point<3> corner2(lx, ly, lz);
    vector<unsigned int> rep;
    rep.push_back(1);
    rep.push_back(2);
    rep.push_back(1);
    dealii::GridGenerator::subdivided_hyper_rectangle(*rect, rep, corner1, corner2, true);
    //// TODO: generate using step_sizes
    return rect;
}

//boost::shared_ptr<Mesh<3> > MixingLayer3D::makeGrid() {
//    //Creation of the principal domain
//    boost::shared_ptr<Mesh<3> > rect = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
//    dealii::Point<3> corner1(-lx/2, -ly/2, -lz/2);
//    dealii::Point<3> corner2(lx/2, ly/2, lz/2);
//    int n_cells_x = 10;
//    int n_cells_y = 10;
//    int n_cells_z = 10;
//    double l_cells_x = lx/n_cells_x;
//    double l_cells_y = ly/n_cells_y;
//    double l_cells_z = lz/n_cells_z;
//    vector<vector<double>> step_sizes;
//    step_sizes.resize(3);
////    for (int ix = 0; ix < n_cells_x; ix++) {
////        step_sizes.at(0).push_back(l_cells_x);
////    }
////    for (int iy = 0; iy < n_cells_y; iy++) {
////        step_sizes.at(1).push_back(l_cells_y);
////    }
////    for (int iz = 0; iz < n_cells_z; iz++) {
////        step_sizes.at(2).push_back(l_cells_z);
////    }
//    step_sizes.at(0) = {lx};
//    step_sizes.at(1) = {ly};
//    step_sizes.at(1) = {ly/2, ly/2};
////    step_sizes.at(1) = {ly/4, ly/8, ly/16, ly/16, ly/16, ly/16, ly/8, ly/4};
//    step_sizes.at(2) = {lz};
//
//    bool colorize = true; 	// set boundary ids automatically to
//    // 0:front; 1:back; 2:bottom; 3:top; 4:left; 5:right;
//    dealii::GridGenerator::subdivided_hyper_rectangle(*rect, step_sizes, corner1, corner2, colorize);
//
//    return rect;
//}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3> > MixingLayer3D::makeBoundaries() {
    // make boundary description
    boost::shared_ptr<BoundaryCollection<3>> boundaries = boost::make_shared<BoundaryCollection<3>>();

    // velocity vector moving forward
    dealii::Vector<double> plusVector(3);
    plusVector[0]=m_U;
    plusVector[1]=0.0;
    plusVector[2]=0.0;

    // velocity vector moving backward
    dealii::Vector<double> minusVector(3);
    minusVector[0]=-m_U;
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