/*
 * MixingLayer3D.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */

#include "MixingLayer3D.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"
#include "natrium/boundaries/DoNothingBoundary.h"
#include "natrium/boundaries/SLFirstOrderBounceBack.h"
#include "natrium/boundaries/ThermalBounceBack.h"
#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include <random>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>
#include <vector>
#include <iostream>
#include "deal.II/grid/grid_in.h"
#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_generator.h"
using namespace std;
using namespace natrium::DealIIExtensions;

double shearlayerthickness = 0.093;

namespace natrium {

MixingLayer3D::MixingLayer3D(double viscosity, size_t refinementLevel, vector<unsigned int> repetitions,
                             double randu_scaling, string randuname, double len_x, double len_y, double len_z,
                             string meshname, double center, double scaling, double dT0, double U, double T,
                             double bc_T, string bc) :
        ProblemDescription<3>(makeGrid(meshname, len_x, len_y, len_z, std::move(repetitions)), viscosity, 1),
                m_initialT(T), m_BCT(bc_T), lx(len_x), ly(len_y), lz(len_z), m_center(center), m_scaling(scaling),
                deltaTheta0(dT0), m_U(U), m_bc(std::move(bc)), m_refinementLevel(refinementLevel) {
    // **** Recommendations for CPU use ****
	/*LOG(BASIC) << "-------------------------------------------------------------" << endl;
	LOG(BASIC) << "**** Recommendations for CPU use ****" << endl;
	double noRepetitions3D = repetitions.at(0) * repetitions.at(1) * repetitions.at(2);
	double noGridPoints = pow( orderOfFiniteElement, 3 ) * pow( 8, refinementLevel +1 ) * noRepetitions3D;
	LOG(BASIC) << "... Computation node details: " << endl;
	LOG(BASIC) << "    - #CPU per node: 12 " << endl;
	LOG(BASIC) << "    - memory per node: 4000 MB " << endl;
	LOG(BASIC) << "... Recommended number of total grid points per node: 10e+6" << endl;
	LOG(BASIC) << "... Recommended number of nodes: " << ceil(noGridPoints/10e+6) << endl;
	LOG(BASIC) << "------------------------------------------------------------" << endl;


	// **** Grid properties ****
	LOG(BASIC) << "**** Grid properties ****" << endl;
	int 	noCellsInXDir	= orderOfFiniteElement * pow( 2, refinementLevel + 1 ) * repetitions.at(0);
	int 	noCellsInYDir	= orderOfFiniteElement * pow( 2, refinementLevel + 1 ) * repetitions.at(1);
	int 	noCellsInZDir	= orderOfFiniteElement * pow( 2, refinementLevel + 1 ) * repetitions.at(2);
	LOG(BASIC) << "... Mesh resolution: " << noCellsInXDir << "x" << noCellsInYDir << "x" << noCellsInZDir << endl;
	LOG(BASIC) << "... Number of total grid points: " << noGridPoints << endl;
    */
    /// apply boundary values
    if (is_MPI_rank_0()) LOG(DETAILED) << "Setting boundaries." << endl;
    setBoundaries(makeBoundaries());
    // apply analytic solution
    if (is_MPI_rank_0()) LOG(DETAILED) << "Setting initial velocity." << endl;
    this->setInitialU(boost::make_shared<InitialVelocity>(this, randu_scaling, randuname, dT0));
    if (is_MPI_rank_0()) LOG(DETAILED) << "Setting initial density." << endl;
    this->setInitialRho(boost::make_shared<InitialDensity>(this));
    if (is_MPI_rank_0()) LOG(DETAILED) << "Setting initial temperature." << endl;
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

MixingLayer3D::~MixingLayer3D() = default;

double MixingLayer3D::InitialVelocity::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component < 3);
    if (m_dT0 == 0) {
        if (component == 0) {
            return 1.;
        } else {
            return 0.;
        }
    }
    double scaling = m_randu_scaling * exp(-pow((x(1))/(2 * m_dT0), 2));
    double rand_u = InterpolateVelocities(x(0), x(1), x(2), component) * scaling;
//    rand_u = 0;
    if (component == 0) {
        return m_flow->m_U * tanh(-x(1)/(2 * m_dT0)) + rand_u;
    } else {
        return rand_u;
    }
}

MixingLayer3D::InitialVelocity::InitialVelocity(natrium::MixingLayer3D *flow, double randu_scaling, string randuname,
                                                double dT0) :
m_flow(flow), m_randu_scaling(randu_scaling), m_dT0(dT0) {
    stringstream filename;
    filename << getenv("NATRIUM_DIR") << "/src/examples/MixingLayer3D/random_u/random_u_" << randuname << ".txt";
    string filestring = filename.str();
    ifstream file(filestring);
    if (is_MPI_rank_0()) LOG(WELCOME) << "Reading random velocity field from " << filestring << endl;
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
                        double tmpk = stod(cell_k);
                        tmpj.push_back(tmpk);
                } tmpi.push_back(tmpj); nz = tmpj.size();
            } tmpdir.push_back(tmpi); ny = tmpi.size();
        } curlOfPsi.push_back(tmpdir); nx = tmpdir.size();
    }
    double lx_u = m_flow->lx;
    double ly_u = m_flow->ly;
    double lz_u = m_flow->lz;
//    double dT0  = m_flow->deltaTheta0;

    if (is_MPI_rank_0()) {
        LOG(DETAILED) << "Creating linspaces x, y, z for interpolation from random_u." << endl
                     << "nx: " << nx << ", ny: " << ny << ", nz: " << nz << endl
                     << "lx: " << lx_u << ", ly: " << ly_u << ", lz: " << lz_u << endl
                     << "lx/dTh0: " << lx_u / m_dT0 << ", ly/dTh0: " << ly_u / m_dT0 << ", lz/dTh0: " << lz_u / m_dT0
                     << endl;
    }

    //// velocity field is scaled to domain
    vector<double> x, y, z;
    double dx, dy, dz;
    double xmin, ymin, zmin;
    xmin = -lx_u / 2;
    ymin = -ly_u / 2;
    zmin = -lz_u / 2;
    dx = lx_u / nx;
    dy = ly_u / ny;
    dz = lz_u / nz;
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
    (void) x;
    return 1.0;// + p / (m_flow->m_cs * m_flow->m_cs);
}
double MixingLayer3D::InitialTemperature::value(const dealii::Point<3>& x, const unsigned int component) const {
    assert(component == 0);
    (void) x;
    return this->m_flow->m_initialT;
}

/**
 * @short create triangulation for Compressible Mixing Layer flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<3> > MixingLayer3D::makeGrid(const string& meshname, double len_x, double len_y, double len_z, std::vector<unsigned int> repetitions) {
    boost::shared_ptr<Mesh<3> > mesh = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);

//    // TODO: generate using step_sizes
//    dealii::Point<3> corner1(-lx/20, -ly/2, -lz/2);
//    dealii::Point<3> corner2(lx/20, ly/2, lz/2);
//    vector<vector<double>> step_sizes(3, vector<double>(20, 10)); // domain is 387 slt, but we need only 20 slt very fine -> divide into 20 areas with decreasing step size
////    step_sizes[0] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
//    step_sizes[1] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.001, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // only central 2 areas are fine
////    step_sizes[2] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
//    dealii::GridGenerator::subdivided_hyper_rectangle(*mesh, step_sizes, corner1, corner2);

    if (meshname == "cube") {
        if (is_MPI_rank_0()) cout << "doing cube with global refinement" << endl;
        boost::shared_ptr<Mesh<3> > cube = boost::make_shared<Mesh<3> >(MPI_COMM_WORLD);
        double lx2, ly2, lz2;
        lx2 = len_x / 2;
        ly2 = len_y / 2;
        lz2 = len_z / 2;
        dealii::Point<3> corner1(-lx2, -ly2, -lz2);
        dealii::Point<3> corner2(lx2, ly2, lz2);
        dealii::GridGenerator::subdivided_hyper_rectangle(*cube, repetitions, corner1, corner2, true);
        return cube;
    } else {
        //Taken from DiamondObstacle2D in step-gridin
        string mesh_filename = "/src/examples/MixingLayer3D/mesh/shearlayer_" + meshname + ".msh";
        if (is_MPI_rank_0()) LOG(WELCOME) << "Reading mesh from " << mesh_filename << endl;
        dealii::GridIn<3> grid_in;
        grid_in.attach_triangulation(*mesh);
        //// Read mesh data from file
        stringstream filename;
        filename << getenv("NATRIUM_DIR") << mesh_filename;
        ifstream file(filename.str().c_str());
        assert(file);
        grid_in.read_msh(file);
    }
    //// calculate boundaries and set boundary ids
    if (is_MPI_rank_0()) LOG(DETAILED) << " dimensions: 3" << endl << " no. of cells: " << mesh->n_active_cells() << endl;
    double minx=0, maxx=0, miny=0, maxy=0, minz=0, maxz=0;
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
    // communicate
    minx = dealii::Utilities::MPI::min_max_avg(minx, MPI_COMM_WORLD).min;
    miny = dealii::Utilities::MPI::min_max_avg(miny, MPI_COMM_WORLD).min;
    minz = dealii::Utilities::MPI::min_max_avg(minz, MPI_COMM_WORLD).min;
    maxx = dealii::Utilities::MPI::min_max_avg(maxx, MPI_COMM_WORLD).max;
    maxy = dealii::Utilities::MPI::min_max_avg(maxy, MPI_COMM_WORLD).max;
    maxz = dealii::Utilities::MPI::min_max_avg(maxz, MPI_COMM_WORLD).max;

    lx = maxx-minx;
    ly = maxy-miny;
    lz = maxz-minz;
    //// set boundary indicators (set top and bottom last to include them in moving bc)
    double rtol = 1e-14;
    for (typename Triangulation<3>::active_cell_iterator cell = mesh->begin_active(); cell != mesh->end(); ++cell) {
        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary()) {
                cell->face(f)->set_all_boundary_ids(42); // just to see after, if any boundaries were unassigned
                Point<3> x = cell->face(f)->center();
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
        LOG(WELCOME) << " domain limits: x in[" << minx << "," << maxx << "], y in [" << miny << "," << maxy << "], z in [" << minz << "," << maxz << "]" << endl;
        LOG(WELCOME) << " boundary indicators: ";
        for (const pair<const types::boundary_id, unsigned int> &pair: boundary_count) {
            LOG(WELCOME) << pair.first << "(" << pair.second << " times) ";
        }
        LOG(WELCOME) << endl;
    }
    mesh->get_boundary_ids();

    return mesh;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<3>> MixingLayer3D::makeBoundaries() {
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

    if (deltaTheta0 == 0) minusVector[0]=m_U;

    // set boundaries on top and bottom to move forward / backward
        if (m_bc == "EQ_BC") {
            boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3>>(2, plusVector, m_BCT));
            boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3>>(3, minusVector, m_BCT));
        }
        if (m_bc == "EQ_DN") {
            boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<3>>(2, plusVector, m_BCT));
            boundaries->addBoundary(boost::make_shared<DoNothingBoundary<3>>(3));
        }
    else if (m_bc == "DN_BC") {
        boundaries->addBoundary(boost::make_shared<DoNothingBoundary<3>>(2));
        boundaries->addBoundary(boost::make_shared<DoNothingBoundary<3>>(3));
    }
    else if (m_bc == "FOBB_BC") {
        boundaries->addBoundary(boost::make_shared<SLFirstOrderBounceBack<3>>(2));
        boundaries->addBoundary(boost::make_shared<SLFirstOrderBounceBack<3>>(3));
    }
    else if (m_bc == "ThBB_BC") {
        boundaries->addBoundary(boost::make_shared<ThermalBounceBack<3>>(2, plusVector, m_initialT));
        boundaries->addBoundary(boost::make_shared<ThermalBounceBack<3>>(3, minusVector, m_initialT));
    }
    else if (m_bc == "VNeq_BC") {
        boundaries->addBoundary(boost::make_shared<VelocityNeqBounceBack<3>>(2, plusVector));
        boundaries->addBoundary(boost::make_shared<VelocityNeqBounceBack<3>>(3, minusVector));
    }
    else if (m_bc == "PP_BC") {
        boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3>>(2, 3, 1, getMesh()));
    }
    if (is_MPI_rank_0()) LOG(DETAILED) << "Boundary condition: " << m_bc << endl;

    // set a boundary between 0 and 1, and 4 and 5, with direction 0 (x) and 2 (z), respectively
    boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3>>(0, 1, 0, getMesh()));
    boundaries->addBoundary(boost::make_shared<PeriodicBoundary<3>>(4, 5, 2, getMesh()));

    // Get the triangulation object (which belongs to the parent class).
    boost::shared_ptr<Mesh<3> > tria_pointer = getMesh();
    return boundaries;
}

} /* namespace natrium */
