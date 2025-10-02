//
// Created by philipp on 9/14/23.
//
#include <random>
#include <vector>
#include "fstream"
#include "stdlib.h"
#include "stdio.h"
#include "iostream"
#include "mpi.h"
#include "natrium/utilities/CommandLineParser.h"
#include </home/philipp/.natrium/libs/fftw/include/fftw3.h>
#define n0 1720
#define n1 172
#define n2 386

using namespace std;

int main(int argc, char** argv) {
    natrium::CommandLineParser parser(argc, argv);
    bool m_print = true;
    auto dirname = parser.getArgument<string>("dirname");


//    const int nx, ny, nz;
//    int kxmax, kymax, kzmax;

    const int nx = 1720; //48; //
    const int ny = 172; //
    const int nz = 386; //48; //
    const int kxmax = nx;
    const int kymax = ny;
    const int kzmax = nz;

    cout << "Recalculating random velocity." << endl;
    if (m_print) cout << "Prints are in " << dirname << "." << endl;

    //// warm up randomness
    static int sseq[mt19937::state_size];
    const static bool once = (srand(std::time(nullptr)), // for true randomness: time(nullptr)
            generate(begin(sseq), end(sseq), rand), true);
    static seed_seq seed_seq(begin(sseq), end(sseq));
    static mt19937 twister(seed_seq);
    // random generator in [-1,1]
    static uniform_real_distribution<double> distr(-1., 1.); // +- Velocity

    vector < vector < vector < vector < double > > > > randomPsi(3);
    for (int dir = 0; dir < 3; dir++) {
        //// set up fftw
        fftw_complex in[nx][ny][nz], out[kxmax][kymax][kzmax];
        fftw_plan p;
        int nx1 = nx, ny1 = ny, nz1 = nz;
        int n[3];
        n[0] = nx; n[1] = ny; n[2] = nz;
        p = fftw_plan_dft(3, *n, in, out, 1, FFTW_FORWARD | FFTW_ESTIMATE);
        p = fftw_plan_dft_3d(n0, n1, n2, in, out, 1, FFTW_FORWARD | FFTW_ESTIMATE);
        p = fftwnd_create_plan(n0, n1, n2, in, out, 1, FFTW_FORWARD | FFTW_ESTIMATE);
        fftwnd_plan pnd = fftw_create_plan(nx1, FFTW_FORWARD, FFTW_ESTIMATE);
//        p = fftw3d_create_plan(nx, ny, nz, FFTW_FORWARD, FFTW_MEASURE);
        fftw_one(p, in, out);
        fftw_destroy_plan(p);

        vector < vector < vector < double > > > tmpdir;
        //// Fill randomPsi with random values
        cout << "Creating random velocity vector potential in direction " << dir << "." << endl;
        for (int xi = 0; xi < nx; xi++) {
            for (int yi = 0; yi < ny; yi++) {
                for (int zi = 0; zi < nz; zi++) {
                    in[xi][yi][zi][0] = distr(twister);
                }
            }
        }
        if (m_print) {
            cout << "Printing random velocity vector potential in direction " << dir << "." << endl;
            ofstream rand_file(dirname + "/random_psi.txt");
            for (int i = 0; i < nx; i++) {
                rand_file << "[";
                for (int j = 0; j < ny; j++) {
                    rand_file << "[";
                    for (int k = 0; k < nz; k++) {
                        if (k > 0) rand_file << ' ';
                        rand_file << tmpdir[i][j][k];
                    } // for all elements in z direction
                    if (j < ny - 1) rand_file << ',';
                } // for all elements in y direction
                if (i < nx - 1) rand_file << ';';
            } // for all elements in x direction
            rand_file << endl;
        }
        cout << "Fourier transforming random velocity vector potential in direction " << dir << "." << endl;
        //// perform dft on randomPsi
        vector < vector < vector < complex < double >> >> psi_hat = Fourier3D(tmpdir);
        if (m_print & is_MPI_rank_0()) {
            cout << "Printing Fourier transform of random velocity vector potential in direction " << dir << "."
                 << endl;
            ofstream dft_file(dirname + "/psi_hat.txt");
            for (int i = 0; i < kxmax; i++) {
                dft_file << "[";
                for (int j = 0; j < kymax; j++) {
                    dft_file << "[";
                    for (int k = 0; k < kzmax; k++) {
                        if (k > 0) dft_file << ' ';
                        dft_file << psi_hat[i][j][k];
                    } // for all elements in z direction
                    if (j < ny - 1) dft_file << ',';
                } // for all elements in y direction
                if (i < nx - 1) dft_file << ';';
            } // for all elements in x direction
            dft_file << endl;
        }
//// multiply in spectral space
        if (is_MPI_rank_0())
            cout << "Scaling Fourier transformed velocity vector potential in direction " << dir << "." << endl;
        double k0 = 23.66 * 0.093; //sqrt(kxmax*kxmax+kymax*kymax+kzmax*kzmax); //
        for (int kx = 0; kx < int(kxmax / 2); kx++) {
            for (int ky = 0; ky < int(kymax / 2); ky++) {
                for (int kz = 0; kz < int(kzmax / 2); kz++) {
                    double k_abs;
                    k_abs = sqrt(kx * kx + ky * ky + kz * kz);
                    psi_hat[kx][ky][kz] *= exp(-2 * k_abs /
                                               k0); // Sarkar: pow(k_abs/k0, 4) * exp(-2 * pow(k_abs / k0, 2)); k0 = 48; // holger: exp(-2 * k_abs / k0); k0 = 23.66*0.093;
                    psi_hat[kxmax - 1 - kx][kymax - 1 - ky][kzmax - 1 - kz] *= exp(-2 * k_abs / k0);
                }
            }
        }
        psi_hat[0][0][0] = 0;
        if (m_print & is_MPI_rank_0()) {
            cout << "Printing scaled Fourier transformed velocity vector potential in direction " << dir << "." << endl;
            ofstream scaling_file(dirname + "/psi_hat_scaled.txt");
            for (int i = 0; i < kxmax; i++) {
                scaling_file << "[";
                for (int j = 0; j < kymax; j++) {
                    scaling_file << "[";
                    for (int k = 0; k < kzmax; k++) {
                        if (k > 0) scaling_file << ' ';
                        scaling_file << psi_hat[i][j][k];
                    } // for all elements in z direction
                    if (j < ny - 1) scaling_file << ',';
                } // for all elements in y direction
                if (i < nx - 1) scaling_file << ';';
            } // for all elements in x direction
            scaling_file << endl;
        }
//// perform inverse dft on psi_hat (directionally)
        if (is_MPI_rank_0())
            cout << "Inversely Fourier transforming velocity vector potential in direction " << dir << "." << endl;
        vector < vector < vector < double>>> psi_i = InverseFourier3D(psi_hat);
        if (m_print & is_MPI_rank_0()) {
            cout << "Printing inverse Fourier transform in direction " << dir << "." << endl;
            ofstream idft_file(dirname + "/psi_idft.txt");
            for (int i = 0; i < nx; i++) {
                idft_file << "[";
                for (int j = 0; j < ny; j++) {
                    idft_file << "[";
                    for (int k = 0; k < nz; k++) {
                        if (k > 0) idft_file << ' ';
                        idft_file << psi_i[i][j][k];
                    } // for all elements in z direction
                    if (j < ny - 1) idft_file << ',';
                } // for all elements in y direction
                if (i < nx - 1) idft_file << ';';
            } // for all elements in x direction
            idft_file << endl;
        }
//// add tmpdir (psix, psiy, or psiz) to randomPsi
        randomPsi.push_back(psi_i);
    } // so: randomPsi = {psix, psiy, psiz} ;

//// calculate gradient using central difference scheme
    if (is_MPI_rank_0()) cout << "Calculating gradients of velocity vector potential." << endl;
    vector < vector < vector < vector < vector < double >>>>> gradient(3,
                                                                       vector < vector < vector < vector <
                                                                       double >> >> (3,
                                                                               vector < vector < vector < double>>>(nx,
            vector < vector < double >> (ny,
                    vector<double>(nz))))); // coordinates are on last three dimensions
    int il, jl, kl, iu, ju, ku;
    for (int dir_psi = 0; dir_psi < 3; dir_psi++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    if (i == 0) { il = 1; } else il = i;
                    if (j == 0) { jl = 1; } else jl = j;
                    if (k == 0) { kl = 1; } else kl = k;
                    if (i == nx - 1) { iu = nx - 2; } else iu = i;
                    if (j == ny - 1) { ju = ny - 2; } else ju = j;
                    if (k == nz - 1) { ku = nz - 2; } else ku = k;
                    gradient[dir_psi][0][i][j][k] =
                            (randomPsi[dir_psi][iu + 1][j][k] - randomPsi[dir_psi][il - 1][j][k]) / (lx / nx);
                    gradient[dir_psi][1][i][j][k] =
                            (randomPsi[dir_psi][i][ju + 1][k] - randomPsi[dir_psi][i][jl - 1][k]) / (ly / ny);
                    gradient[dir_psi][2][i][j][k] =
                            (randomPsi[dir_psi][i][j][ku + 1] - randomPsi[dir_psi][i][j][kl - 1]) / (lz / nz);
                }
            }
        }
    } // so: gradient = {gradient_psix, gradient_psiy, gradient_psiz} and gradient_psin = { dpsin/dx, dpsin/dy, dpsin/dz }
    if (m_print & is_MPI_rank_0()) {
        cout << "Printing gradients of velocity vector potential." << endl;
        ofstream grad_file(dirname + "/psi_grad.txt");
        for (int dir_psi = 0; dir_psi < 3; dir_psi++) {
            grad_file << "direction: " << to_string(dir_psi) << endl << "[";
            for (int dir_grad = 0; dir_grad < 3; dir_grad++) {
                grad_file << "[";
                for (int i = 0; i < nx; i++) {
                    grad_file << "[";
                    for (int j = 0; j < ny; j++) {
                        grad_file << "[";
                        for (int k = 0; k < nz; k++) {
                            if (k > 0) grad_file << ' ';
                            grad_file << gradient[dir_psi][dir_grad][i][j][k];
                        } // for all elements in z direction
                        if (j < ny - 1) grad_file << ',';
                    } // for all elements in y direction
                    if (i < nx - 1) grad_file << ';';
                } // for all elements in x direction
                if (dir_grad < 2) grad_file << "]";
            }
            grad_file << endl;
        }
    }
//// calculate curl using gradient values
    if (is_MPI_rank_0()) cout << "Calculating curl of velocity vector potential." << endl;
    int m, n; // for indices of cross-product
    for (int dir_curl = 0; dir_curl < 3; dir_curl++) {
        vector < vector < vector < double > > > tmpdir;
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
            vector <vector<double>> tmpi;
            for (int j = 0; j < ny; j++) {
                vector<double> tmpj;
                for (int k = 0; k < nz; k++) {
                    double tmpk;
                    tmpk = (gradient[m][n][i][j][k] - gradient[n][m][i][j][k]) / 4;
                    tmpj.push_back(tmpk);
                }
                tmpi.push_back(tmpj);
            }
            tmpdir.push_back(tmpi);
        }
        curlOfPsi.push_back(tmpdir);
    } // so: curlOfPsi = {ux, uy, uz}
    if (m_print & is_MPI_rank_0()) {
        cout << "Printing velocity vector matrix." << endl;
        ofstream u_file(dirname + "/random_u.txt");
        for (int dir_psi = 0; dir_psi < 3; dir_psi++) {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        if (k > 0) u_file << ' ';
                        u_file << curlOfPsi[dir_psi][i][j][k];
                    } // for all elements in z direction
                    if (j < ny - 1) u_file << ',';
                } // for all elements in y direction
                if (i < nx - 1) u_file << ';';
            } // for all elements in x direction
            u_file << endl;
        }
    }
}
