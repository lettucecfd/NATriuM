/**
 * @file DiamondObstacle2D.cpp
 * @short Flow around a diamond-shaped obstacle
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DiamondObstacle2D.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/grid/grid_in.h"

#include "natrium/boundaries/VelocityNeqBounceBack.h"
#include "natrium/boundaries/DoNothingBoundary.h"
#include "natrium/boundaries/SLEquilibriumBoundary.h"

#include "natrium/utilities/Math.h"

namespace natrium {



DiamondObstacle2D::DiamondObstacle2D(double velocity, double viscosity, size_t refinementLevel, int aoa, string foilname):
		ProblemDescription<2>(makeGrid(refinementLevel, aoa, foilname), viscosity, 1.0), m_meanInflowVelocity(velocity),
        m_refinementLevel(refinementLevel) {

	/// apply boundary values
	setBoundaries(makeBoundaries());
//	double Fx = 8 * velocity * viscosity;
//  if(is_MPI_rank_0()) LOG(DETAILED) << "F: " << Fx << endl;
//	dealii::Tensor<1, 2> F;
//	F[0] = Fx;
//  F[0] = Fx*0.04;
//  setExternalForce(boost::make_shared<ConstantExternalForce<2> >(F));
    this->setInitialU(boost::make_shared<InitialVelocity>(this));
    this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

DiamondObstacle2D::~DiamondObstacle2D() {}

/**
 * @short Read the mesh given in deal.ii's step-35
 * @return shared pointer to a triangulation instance
 */
double DiamondObstacle2D::InitialVelocity::value(const dealii::Point<2>& x,
		const unsigned int component) const {
    assert (component < 2);
    if (component == 0) return this->m_flow->m_meanInflowVelocity;
	else return this->m_flow->m_meanInflowVelocity*sin(x[1]*1.7+0.1)*0.1;
}

double DiamondObstacle2D::InitialDensity::value(const dealii::Point<2>& x, const unsigned int component) const {
    (void) x;
    (void) component;
    return 1.0;
}

double DiamondObstacle2D::InitialTemperature::value(const dealii::Point<2>& x, const unsigned int component) const {
    (void) x;
    (void) component;
    return 1.0;
}

boost::shared_ptr<Mesh<2> > DiamondObstacle2D::makeGrid(size_t refinementLevel, int aoa, string foilname) {
    (void) refinementLevel;
	//Read in grid
	//Taken from step-35 in deal.ii
    if (is_MPI_rank_0()) LOG(WELCOME) << "Mesh name set to " << foilname << endl;
	dealii::GridIn<2> grid_in;
	boost::shared_ptr<Mesh<2>> mesh = boost::make_shared<Mesh<2>>(MPI_COMM_WORLD);
	grid_in.attach_triangulation(*mesh);
	{
		std::stringstream filename;
		filename << getenv("NATRIUM_DIR") << "/src/examples/step-grid-in/mesh/" << foilname << "/NACA0012_" << aoa << "deg.msh"; // "/src/examples/step-grid-in/Archive/naca0012_supersonic.msh";//
		std::ifstream file(filename.str().c_str());
        if (is_MPI_rank_0()) LOG(WELCOME) << "Reading mesh from " << filename.str() << endl;
		assert(file);
        grid_in.read_msh(file);
	}
	// mesh->refine_global (refinementLevel);
    if (is_MPI_rank_0()) LOG(WELCOME) << "Getting boundary IDs." << endl;
    mesh->get_boundary_ids();
	return mesh;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > DiamondObstacle2D::makeBoundaries() {
    if (is_MPI_rank_0()) LOG(WELCOME) << "Making boundaries." << endl;
	// make boundary description
	boost::shared_ptr<BoundaryCollection<2>> boundaries = boost::make_shared<BoundaryCollection<2>>();
	dealii::Vector<double> inflow(2);
	inflow[0]=m_meanInflowVelocity;
	inflow[1]=0.0;

    dealii::Tensor<1, 2> u;
    u[0] = inflow(0);
    u[1] = inflow(1);

    dealii::Tensor<1, 2> profile;
        profile[0] = 0.0;
        profile[1] = 0.0;

	std::vector<double> oneVector = {m_meanInflowVelocity, 1.0};
//    oneVector[0]=m_meanInflowVelocity;
//    oneVector[1]=1.0;
	boost::shared_ptr<dealii::Function<2>> boundary_density = boost::make_shared<dealii::Functions::ConstantFunction<2>> (1.0);
	boost::shared_ptr<dealii::Function<2>> boundary_velocity = boost::make_shared<InflowVelocity> (m_meanInflowVelocity);

    boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<2>>(300, inflow));
//    boundaries->addBoundary(boost::make_shared<SLEquilibriumBoundary<2>>(302, inflow)); // outflow
	boundaries->addBoundary(boost::make_shared<DoNothingBoundary<2>>(302)); // outflow
	boundaries->addBoundary(boost::make_shared<DoNothingBoundary<2>>(301)); // sponge
    boundaries->addBoundary(boost::make_shared<VelocityNeqBounceBack<2>>(303, profile));

	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
