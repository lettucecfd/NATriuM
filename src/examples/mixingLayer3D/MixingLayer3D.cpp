/*
 * MixingLayer3D.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */

#include "MixingLayer3D.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"
#include <random>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>
#include "deal.II/grid/grid_in.h"
#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_tools.h"
using namespace std;
using namespace natrium::DealIIExtensions;

double shearlayerthickness = 0.093;
//double k0 = 23.66 * shearlayerthickness; // peak wave number
//int n = 5;
//int kmax = pow(2, n); // [1] C. Pantano and S. Sarkar, “A study of compressibility effects in the high-speed turbulent shear layer using direct simulation,” J. Fluid Mech., vol. 451, pp. 329–371, Jan. 2002, doi: 10.1017/S0022112001006978.
// kmax = 32
//int npoints = 32; // number of points in shortest axis of velocity field (lz, presumably)
//bool print = false;
//bool recalculate_psi = true;

namespace natrium {

MixingLayer3D::MixingLayer3D(double viscosity, size_t refinementLevel, bool squash, bool print, bool recalculate,
                             string dirName, double U) :
ProblemDescription<3>(makeGrid(), viscosity, 1), m_squash(squash), m_U(U),
lx(1720*shearlayerthickness), ly(387*shearlayerthickness), lz(172*shearlayerthickness), m_refinementLevel(refinementLevel) {
//    if (m_refinementLevel > 4) { print = false; }
    /// apply boundary values
    setBoundaries(makeBoundaries());
    // apply analytic solution
    this->setInitialU(boost::make_shared<InitialVelocity>(this, print, recalculate, dirName));
    this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

MixingLayer3D::~MixingLayer3D() = default;

double MixingLayer3D::InitialVelocity::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component < 3);
    double scaling = 0.01 * exp(-pow((x(1))/(2 * shearlayerthickness), 2));
    double rand_u = InterpolateVelocities(x(0), x(1), x(2), component) * scaling;
//    rand_u = 0;
    if (component == 0) {
        return tanh(-x(1)/(2 * shearlayerthickness)) + rand_u;
    } else {
        return rand_u;
    }
}

MixingLayer3D::InitialVelocity::InitialVelocity(natrium::MixingLayer3D *flow, bool print, bool recalculate, string dirName) :
m_flow(flow), lx(flow->lx), ly(flow->ly), lz(flow->lz), m_print(print), m_recalculate(recalculate) {
//    int kmax = 48;//pow(2, flow->m_refinementLevel);
//    nx = 1720/2; // pow(2, 5);
//    ny = 387/2; //ny;
//    nz = 172/2; //nz;
//    kxmax = nx / 5;
//    kymax = ny / 5;
//    kzmax = nz / 5;

    if (m_recalculate) {
        nx = pow(2, 6); //48; //
        ny = pow(2, 6); //48; //
        nz = pow(2, 6); //48; //
        kxmax = nx;
        kymax = ny;
        kzmax = nz;

        if (is_MPI_rank_0()) {
            cout << "Recalculating random velocity." << endl;
            if (m_print) cout << "Prints are in " << dirName << "." << endl;
        }

        //// warm up randomness
        static int sseq[ mt19937::state_size ] ;
        const static bool once = ( srand( std::time(nullptr)), // for true randomness: time(nullptr)
                generate( begin(sseq), end(sseq), rand ), true ) ;
        static seed_seq seed_seq( begin(sseq), end(sseq) ) ;
        static mt19937 twister(seed_seq) ;
        // random generator in [-1,1]
        static uniform_real_distribution<double> distr(-1., 1.) ; // +- Velocity

        randomPsi.reserve(3);
        for (int dir = 0; dir < 3; dir++) { vector< vector< vector<double> > > tmpdir;
            //// Fill randomPsi with random values
            if (is_MPI_rank_0()) cout << "Creating random velocity vector potential in direction " << dir << "." << endl;
            for (int xi = 0; xi < nx; xi++) { vector<vector<double> > tmpi;
                for (int yi = 0; yi < ny; yi++) { vector<double> tmpj;
                    for (int zi = 0; zi < nz; zi++) { double tmpk;
                        tmpk = distr(twister);
                        tmpj.push_back(tmpk);
                    } tmpi.push_back(tmpj);
                } tmpdir.push_back(tmpi);
            }
            if (m_print & is_MPI_rank_0()) {
                cout << "Printing random velocity vector potential in direction " << dir << "." << endl;
                ofstream rand_file(dirName + "/random_psi.txt");
                for (int i = 0; i < nx; i++) { rand_file << "[";
                    for (int j = 0; j < ny; j++) { rand_file << "[";
                        for (int k = 0; k < nz; k++) {
                            if (k > 0) rand_file << ' ';
                            rand_file << tmpdir[i][j][k];
                        } // for all elements in z direction
                        if (j<ny-1) rand_file << ',';
                    } // for all elements in y direction
                    if (i<nx-1) rand_file << ';';
                } // for all elements in x direction
                rand_file << endl;
            }
            if (is_MPI_rank_0()) cout << "Fourier transforming random velocity vector potential in direction " << dir << "." << endl;
            //// perform dft on randomPsi
            vector< vector< vector<complex<double>>>> psi_hat = Fourier3D(tmpdir);
            if (m_print & is_MPI_rank_0()) {
                cout << "Printing Fourier transform of random velocity vector potential in direction " << dir << "." << endl;
                ofstream dft_file(dirName + "/psi_hat.txt");
                for (int i = 0; i < kxmax; i++) { dft_file << "[";
                    for (int j = 0; j < kymax; j++) { dft_file << "[";
                        for (int k = 0; k < kzmax; k++) {
                            if (k > 0) dft_file << ' ';
                            dft_file << psi_hat[i][j][k];
                        } // for all elements in z direction
                        if (j<ny-1) dft_file << ',';
                    } // for all elements in y direction
                    if (i<nx-1) dft_file << ';';
                } // for all elements in x direction
                dft_file << endl;
            }
            //// multiply in spectral space
            if (is_MPI_rank_0()) cout << "Scaling Fourier transformed velocity vector potential in direction " << dir << "." << endl;
            double k0 = 23.66*0.093; //sqrt(kxmax*kxmax+kymax*kymax+kzmax*kzmax); //
            for (int kx = 0; kx < int(kxmax/2); kx++) {
                for (int ky = 0; ky < int(kymax/2); ky++) {
                    for (int kz = 0; kz < int(kzmax/2); kz++) {
                        double k_abs;
                        k_abs = sqrt(kx * kx + ky * ky + kz * kz);
                        psi_hat[kx][ky][kz] *= exp(-2 * k_abs / k0); // Sarkar: pow(k_abs/k0, 4) * exp(-2 * pow(k_abs / k0, 2)); k0 = 48; // holger: exp(-2 * k_abs / k0); k0 = 23.66*0.093;
                        psi_hat[kxmax-1 - kx][kymax-1 - ky][kzmax-1 - kz] *= exp(-2 * k_abs / k0);
                    } } }
            psi_hat[0][0][0] = 0;
            if (m_print & is_MPI_rank_0()) {
                cout << "Printing scaled Fourier transformed velocity vector potential in direction " << dir << "." << endl;
                ofstream scaling_file(dirName + "/psi_hat_scaled.txt");
                for (int i = 0; i < kxmax; i++) { scaling_file << "[";
                    for (int j = 0; j < kymax; j++) { scaling_file << "[";
                        for (int k = 0; k < kzmax; k++) {
                            if (k > 0) scaling_file << ' ';
                            scaling_file << psi_hat[i][j][k];
                        } // for all elements in z direction
                        if (j<ny-1) scaling_file << ',';
                    } // for all elements in y direction
                    if (i<nx-1) scaling_file << ';';
                } // for all elements in x direction
                scaling_file << endl;
            }
            //// perform inverse dft on psi_hat (directionally)
            if (is_MPI_rank_0()) cout << "Inversely Fourier transforming velocity vector potential in direction " << dir << "." << endl;
            vector<vector<vector<double>>> psi_i = InverseFourier3D(psi_hat);
            if (m_print & is_MPI_rank_0()) {
                cout << "Printing inverse Fourier transform in direction " << dir << "." << endl;
                ofstream idft_file(dirName + "/psi_idft.txt");
                for (int i = 0; i < nx; i++) { idft_file << "[";
                    for (int j = 0; j < ny; j++) { idft_file << "[";
                        for (int k = 0; k < nz; k++) {
                            if (k > 0) idft_file << ' ';
                            idft_file << psi_i[i][j][k];
                        } // for all elements in z direction
                        if (j<ny-1) idft_file << ',';
                    } // for all elements in y direction
                    if (i<nx-1) idft_file << ';';
                } // for all elements in x direction
                idft_file << endl;
            }
            //// add tmpdir (psix, psiy, or psiz) to randomPsi
            randomPsi.push_back(psi_i);
        } // so: randomPsi = {psix, psiy, psiz} ;

        //// calculate gradient using central difference scheme
        if (is_MPI_rank_0()) cout << "Calculating gradients of velocity vector potential." << endl;
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
                        gradient[dir_psi][0][i][j][k] = (randomPsi[dir_psi][iu+1][j][k] - randomPsi[dir_psi][il-1][j][k]) / (lx/nx);
                        gradient[dir_psi][1][i][j][k] = (randomPsi[dir_psi][i][ju+1][k] - randomPsi[dir_psi][i][jl-1][k]) / (ly/ny);
                        gradient[dir_psi][2][i][j][k] = (randomPsi[dir_psi][i][j][ku+1] - randomPsi[dir_psi][i][j][kl-1]) / (lz/nz);
                    }}}} // so: gradient = {gradient_psix, gradient_psiy, gradient_psiz} and gradient_psin = { dpsin/dx, dpsin/dy, dpsin/dz }
        if (m_print & is_MPI_rank_0()) {
            cout << "Printing gradients of velocity vector potential." << endl;
            ofstream grad_file(dirName + "/psi_grad.txt");
            for (int dir_psi = 0; dir_psi < 3; dir_psi++) { grad_file << "direction: " << to_string(dir_psi) << endl << "[";
                for (int dir_grad = 0; dir_grad < 3; dir_grad++) { grad_file << "[";
                    for (int i = 0; i < nx; i++) { grad_file << "[";
                        for (int j = 0; j < ny; j++) { grad_file << "[";
                            for (int k = 0; k < nz; k++) {
                                if (k > 0) grad_file << ' ';
                                grad_file << gradient[dir_psi][dir_grad][i][j][k];
                            } // for all elements in z direction
                            if (j<ny-1) grad_file << ',';
                        } // for all elements in y direction
                        if (i<nx-1) grad_file << ';';
                    } // for all elements in x direction
                    if (dir_grad<2) grad_file << "]";
                } grad_file << endl;
            }
        }
        //// calculate curl using gradient values
        if (is_MPI_rank_0()) cout << "Calculating curl of velocity vector potential." << endl;
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
        if (m_print & is_MPI_rank_0()) {
            cout << "Printing velocity vector matrix." << endl;
            ofstream u_file(dirName + "/random_u.txt");
            for (int dir_psi = 0; dir_psi < 3; dir_psi++) {
                for (int i = 0; i < nx; i++) {
                    for (int j = 0; j < ny; j++) {
                        for (int k = 0; k < nz; k++) {
                            if (k > 0) u_file << ' ';
                            u_file << curlOfPsi[dir_psi][i][j][k];
                        } // for all elements in z direction
                        if (j<ny-1) u_file << ',';
                    } // for all elements in y direction
                    if (i<nx-1) u_file << ';';
                } // for all elements in x direction
                u_file << endl;
            }
        }
    }
    else {
        stringstream filename;
        filename << getenv("NATRIUM_DIR") << "/src/examples/mixingLayer3D/random_u_kmax64.txt";
        string filestring = filename.str();
        ifstream file(filestring);
        if (is_MPI_rank_0()) cout << "Reading initial velocities from " << filestring << endl;
        string line;
        while (getline(file, line)) {
            stringstream linestream(line);
            string cell_i;
            vector<vector<vector<double>>> tmpdir;
            while(getline(linestream, cell_i, ';')) {
                stringstream cell_i_stream(cell_i);
                string cell_j;
                vector<vector<double>> tmpi;
                while(getline(cell_i_stream,cell_j,',')) {
                    stringstream cell_j_stream(cell_j);
                    string cell_k;
                    vector<double> tmpj;
                    while(getline(cell_j_stream, cell_k, ' ')) {
//                        if (!cell_k.empty()) {
                            double tmpk = stod(cell_k);
                            tmpj.push_back(tmpk);
//                        }
                    } tmpi.push_back(tmpj); nz = tmpj.size();
                } tmpdir.push_back(tmpi); ny = tmpi.size();
            } curlOfPsi.push_back(tmpdir); nx = tmpdir.size();
        }
        if (is_MPI_rank_0()) cout << "nx: " << nx << ", ny: " << ny << ", nz: " << nz << endl;
    }

    if (is_MPI_rank_0()) cout << "Creating linspaces x, y, z." << endl;
    vector<double> x, y, z;
    double dx, dy, dz;
    double xmin, ymin, zmin;
    xmin = -lx / 2;
    ymin = -ly / 2;
    zmin = -lz / 2;
    dx = lx / nx;
    dy = ly / ny;
    dz = lz / nz;
    double linvalue;
    linvalue = xmin; for (int i = 0; i < nx; i++) { x.push_back(linvalue); linvalue += dx; }
    linvalue = ymin; for (int i = 0; i < ny; i++) { y.push_back(linvalue); linvalue += dy; }
    linvalue = zmin; for (int i = 0; i < nz; i++) { z.push_back(linvalue); linvalue += dz; }
    xvec.insert(xvec.end(), x.begin(), x.end());
    yvec.insert(yvec.end(), y.begin(), y.end());
    zvec.insert(zvec.end(), z.begin(), z.end());

    auto xbounds = minmax_element(x.begin(), x.end());
    minx = *xbounds.first; maxx = *xbounds.second;
    auto ybounds = minmax_element(y.begin(), y.end());
    miny = *ybounds.first; maxy = *ybounds.second;
    auto zbounds = minmax_element(z.begin(), z.end());
    minz = *zbounds.first; maxz = *zbounds.second;
}

vector<vector<vector<complex<double>>>>
MixingLayer3D::InitialVelocity::Fourier3D(const vector<vector<vector<double>>> &in) const {
    vector<vector<vector<complex<double>>>> out(
            kxmax, vector<vector<complex<double>>>(
                    kymax, vector<complex<double>>(
                            kzmax, {0.0, 0.0})));
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

vector<vector<vector<double>>>
MixingLayer3D::InitialVelocity::InverseFourier3D(const vector<vector<vector<complex<double>>>> &in) const {
    vector<vector<vector<double>>> out(nx,
        vector<vector<double>>(ny,
            vector<double>(nz, 0.0)));
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
    // set coordinates to boundaries, if point is outside
    xq = max(minx, min(xq, maxx));
    yq = max(miny, min(yq, maxy));
    zq = max(minz, min(zq, maxz));
    // find upper and lower bound for interpolation
    auto xupper = upper_bound(xvec.cbegin(), xvec.cend(), xq);
    int x1 = (xupper == xvec.cend()) ? xupper - xvec.cbegin() - 1 : xupper - xvec.cbegin();
    int x0 = x1 - 1;
    auto yupper = upper_bound(yvec.cbegin(), yvec.cend(), yq);
    int y1 = (yupper == yvec.cend()) ? yupper - yvec.cbegin() - 1 : yupper - yvec.cbegin();
    auto y0 = y1 - 1;
    auto zupper = upper_bound(zvec.cbegin(), zvec.cend(), zq);
    int z1 = (zupper == zvec.cend()) ? zupper - zvec.cbegin() - 1 : zupper - zvec.cbegin();
    auto z0 = z1 - 1;
    // get relative distances for weighting
    double xd = (xq - xvec[x0])/(xvec[x1] - xvec[x0]);
    double yd = (yq - yvec[y0])/(yvec[y1] - yvec[y0]);
    double zd = (zq - zvec[z0])/(zvec[z1] - zvec[z0]);
    // get surrounding values
    double c000 = curlOfPsi[dim][x0][y0][z0];
    double c010 = curlOfPsi[dim][x0][y1][z0];
    double c100 = curlOfPsi[dim][x1][y0][z0];
    double c110 = curlOfPsi[dim][x1][y1][z0];
    double c001 = curlOfPsi[dim][x0][y0][z1];
    double c011 = curlOfPsi[dim][x0][y1][z1];
    double c101 = curlOfPsi[dim][x1][y0][z1];
    double c111 = curlOfPsi[dim][x1][y1][z1];
    // calculate interpolations stepwise (x, then y, then z)
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

    //Taken from DiamondObstacle2D in step-gridin
    dealii::GridIn<3> grid_in;
    boost::shared_ptr<Mesh<3> > mesh = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
    grid_in.attach_triangulation(*mesh);

    //// Read mesh data from file
    stringstream filename;
    filename << getenv("NATRIUM_DIR") << "/src/examples/mixingLayer3D/shearlayer_wideCentre.msh";
    ifstream file(filename.str().c_str());
    assert(file);
    grid_in.read_msh(file);
    if (is_MPI_rank_0()) cout << "Imported mesh info:" << endl << " dimension: 3" << endl << " no. of cells: " << mesh->n_active_cells() << endl;

    double minx, maxx, miny, maxy, minz, maxz;
    //// get minimum and maximum coordinates
    for (typename Triangulation<3>::active_cell_iterator cell = mesh->begin_active(); cell != mesh->end(); ++cell) {
        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary()) {
                Point<3> x = cell->face(f)->center();
                minx = min(minx, x[0]);
                maxx = max(maxx, x[0]);
                miny = min(miny, x[1]);
                maxy = max(maxy, x[1]);
                minz = min(minz, x[2]);
                maxz = max(maxz, x[2]);
            }
        }
    }
    //// set boundary indicators (set top and bottom last to include them in moving bc)
    double rtol = 1e-14;
    for (typename Triangulation<3>::active_cell_iterator cell = mesh->begin_active(); cell != mesh->end(); ++cell) {
        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary()) {
                cell->face(f)->set_all_boundary_ids(42); // just to see after, if any boundaries were unassigned
                Point<3> x = cell->face(f)->center();
//                if (x[0] == maxx) // front (x)
//                    cell->face(f)->set_all_boundary_ids(0);
//                if (x[0] == minx) // back (x)
//                    cell->face(f)->set_all_boundary_ids(1);
//                if (x[2] == maxz) // left (z)
//                    cell->face(f)->set_all_boundary_ids(4);
//                if (x[2] == minz) // right (z)
//                    cell->face(f)->set_all_boundary_ids(5);
//                if (x[1] == maxy) // top (y)
//                    cell->face(f)->set_all_boundary_ids(2);
//                if (x[1] == miny) // bottom (y)
//                    cell->face(f)->set_all_boundary_ids(3);
                if (abs(x[0] - maxx)/abs(maxx) < rtol) // front (x)
                    cell->face(f)->set_all_boundary_ids(0);
                if (abs(x[0] - minx)/abs(minx) < rtol) // back (x)
                    cell->face(f)->set_all_boundary_ids(1);
                if (abs(x[2] - maxz)/abs(maxz) < rtol) // left (z)
                    cell->face(f)->set_all_boundary_ids(4);
                if (abs(x[2] - minz)/abs(minz) < rtol) // right (z)
                    cell->face(f)->set_all_boundary_ids(5);
                if (abs(x[1] - miny)/abs(miny) < rtol) // bottom (y)
                    cell->face(f)->set_all_boundary_ids(2);
                if (abs(x[1] - maxy)/abs(maxy) < rtol) // top (y)
                    cell->face(f)->set_all_boundary_ids(3);
            }
        }
    }

    if (is_MPI_rank_0()) {
        //// Print boundary indicators
        map<types::boundary_id, unsigned int> boundary_count;
        for (auto cell: mesh->active_cell_iterators()) {
            for (unsigned int face = 0; face < GeometryInfo<3>::faces_per_cell; ++face) {
                if (cell->face(face)->at_boundary())
                    boundary_count[cell->face(face)->boundary_id()]++;
            }
        }
        cout << " boundary indicators: ";
        for (const pair<const types::boundary_id, unsigned int> &pair: boundary_count) {
            cout << pair.first << "(" << pair.second << " times) ";
        }
        cout << endl;
    }
    mesh->get_boundary_ids();
    return mesh;
}

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

//boost::shared_ptr<Mesh<3> > MixingLayer3D::makeGrid() {
//    //Creation of the principal domain
//    boost::shared_ptr<Mesh<3> > mesh = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
//    double lx = 1720 * shearlayerthickness / 2; // 1720*0.093/2 +-80(160) // dx 2.24 -> nx 70; 64; 2^6 //
//    double ly = 387 * shearlayerthickness / 2; // 387*0.093/2 +-18(36)// dy 1.73..5.32 -> ny 21..7; 16; 2^4 //
//    double lz = 172 * shearlayerthickness / 2; // 172*0.093/2 +-8(16) // dz 1.08..2.69 -> nz 14.8..6; 16; 2^4 //
//    dealii::Point<3> corner1(-lx, -ly, -lz);
//    dealii::Point<3> corner2(lx, ly, lz);
//    vector<unsigned int> rep;
//    rep.push_back(1);
//    rep.push_back(2);
//    rep.push_back(1);
//    dealii::GridGenerator::subdivided_hyper_rectangle(*mesh, rep, corner1, corner2, true);

//    //// add periodicity
//    vector<GridTools::PeriodicFacePair<Mesh<3>::cell_iterator >> periodicity_vector;
//    //    Tensor<1, 3> offset = dealii::Tensor<1, 3 >();
//    //    FullMatrix<double> matrix;
//    const types::boundary_id b_id1=0, b_id2=1;
//    const unsigned int direction = 1;
//    //    GridTools::collect_periodic_faces(*mesh, b_id1, b_id2, direction, periodicity_vector, offset, matrix);
//    GridTools::collect_periodic_faces(*mesh, b_id1, b_id2, direction, periodicity_vector);
//    for (const auto & pair : periodicity_vector) {
//        cout << pair.face_idx << endl;
//    }
//    mesh->add_periodicity(periodicity_vector);

// mesh->refine_global (refinementLevel);

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