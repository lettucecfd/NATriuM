/**
 * @file ShockVortexInteraction.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "ShockVortexInteraction.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "natrium/boundaries/PeriodicBoundary.h"
#include "natrium/boundaries/VelocityNeqBounceBack.h"

#include "natrium/utilities/Math.h"

namespace natrium {

ShockVortexInteraction::ShockVortexInteraction(double viscosity, size_t refinement_level, double u0,
		double kappa, double Ma_v, double perturbation, double trafo_x, double trafo_y) :
		ProblemDescription<2>(makeGrid(), viscosity, 1.0), m_u0(u0), m_kappa(
				kappa), m_refinementLevel(refinement_level), m_perturbation(perturbation),
       	              		m_trafoX(trafo_x), m_trafoY(trafo_y), m_Ma_v(Ma_v)	{
	assert(trafo_x >=0);
	assert(trafo_x < 1);
	assert(trafo_y >=0);
	assert(trafo_y < 1);
	/// apply boundary values
	setBoundaries(makeBoundaries());
	// apply initial and analytical solution
	this->setInitialU(boost::make_shared<InitialVelocity>(this));
	this->setInitialRho(boost::make_shared<InitialDensity>(this));
    this->setInitialT(boost::make_shared<InitialTemperature>(this));
}

ShockVortexInteraction::~ShockVortexInteraction() {
}

    double ShockVortexInteraction::InitialVelocity::value(const dealii::Point<2> &x,
                                                          const unsigned int component) const {
        double u_v = m_flow->m_Ma_v / sqrt(3.0) * sqrt(1.4);
        double offset = m_flow->m_shockPosition; // location of the shock
        double left = -m_flow->m_machNumberLeft / sqrt(3.0) * sqrt(1.4 * m_flow->m_temperatureLeft);
        double right = -m_flow->m_machNumberRight / sqrt(3.0) * sqrt(1.4);
        double steepness = m_flow->m_shockSteepness;

        std::pair<double,double> coeff = m_flow->calcOffsets(left, right);


        double return_value = 0.0;
        if (component == 0) {
            return_value = (tanh(steepness * (x[0] - offset)) + coeff.first) * coeff.second;
        }
        if (component == 1) {
            return_value = 0.0;
        }


        double x_rel = x(0) - (m_flow->m_shockPosition + m_flow->m_vortexOffset);
        double y_rel = x(1) - 12.0;
        double r = sqrt(x_rel * x_rel + y_rel * y_rel);
        double sinalpha = y_rel / r;
        double cosalpha = x_rel / r;

        if (r <= 4.0) {
            if (component == 0) {
                return_value -= r * u_v * sinalpha * exp((1.0 - r * r) / 2.0);
            }

            if (component == 1) {
                return_value += r * u_v * cosalpha * exp((1.0 - r * r) / 2.0);
            }

        }

        return return_value;

    }


    double
    ShockVortexInteraction::InitialDensity::value(const dealii::Point<2> &x, const unsigned int component) const {

        double x_rel = x(0) - (m_flow->m_shockPosition + m_flow->m_vortexOffset);
        double y_rel = x(1) - 12.0;
        double r = sqrt(x_rel * x_rel + y_rel * y_rel);
        double offset = m_flow->m_shockPosition; // location of the shock
        double steepness = m_flow->m_shockSteepness;
        std::pair<double,double> coeff = m_flow->calcOffsets(left, right);

        double return_value = (tanh(steepness * (x[0] - offset)) + tanh_c) * tanh_d;

        double Ma_v = m_flow->m_Ma_v;
        if (r <= 4.0) {
            return_value = pow(1. - ((1.4 - 1.) / 2. * Ma_v * Ma_v * exp(1.0 - r * r)), (1. / (1.4 - 1.)));
        }
        return return_value;
    }




double ShockVortexInteraction::InitialTemperature::value(const dealii::Point<2>& x, const unsigned int component) const {

    double offset = m_flow->m_shockPosition; // location of the shock
    double steepness = m_flow->m_shockSteepness;
    auto[tanh_c, tanh_d] = m_flow->calcOffsets(m_flow->m_temperatureLeft, m_flow->m_temperatureRight);
    std::pair<double,double> coeff = m_flow->calcOffsets(left, right);

    double Ma_v = m_flow->m_Ma_v;
    double x_rel = x(0)  - (m_flow->m_shockPosition+m_flow->m_vortexOffset);
    double y_rel = x(1) - 12.0;
    double r = sqrt(x_rel*x_rel+y_rel*y_rel);
    if(r<=4.0)
    {
       return_value = (1.-(1.4-1.)/2.*Ma_v*Ma_v*exp(1.0-r*r));
    }


return return_value;
}


/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
boost::shared_ptr<Mesh<2> > ShockVortexInteraction::makeGrid() {
	//Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >(
	MPI_COMM_WORLD);
#else
	boost::shared_ptr<Mesh<2> > rect = boost::make_shared<Mesh<2> >();
#endif
	const dealii::Point<2> left = {0.0,0.0};
        const dealii::Point<2> right = {72.0,24.0};
    const std::vector <unsigned int>& reps = {1,1};


	dealii::GridGenerator::subdivided_hyper_rectangle(*rect, reps, left, right, true);
	//dealii::GridGenerator::hyper_cube(*rect, 0, 1);
	// Assign boundary indicators to the faces of the "parent cell"
	Mesh<2>::active_cell_iterator cell = rect->begin_active();
	cell->face(0)->set_all_boundary_ids(0);  // left
	cell->face(1)->set_all_boundary_ids(1);  // right
	cell->face(2)->set_all_boundary_ids(2);  // top
	cell->face(3)->set_all_boundary_ids(3);  // bottom

	return rect;
}

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
boost::shared_ptr<BoundaryCollection<2> > ShockVortexInteraction::makeBoundaries() {

    // make boundary description
    boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
                    BoundaryCollection<2> >();
    numeric_vector zeroVelocity(2);

    boundaries->addBoundary(
                    boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
    boundaries->addBoundary(
                    boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));


	// Get the triangulation object (which belongs to the parent class).
	boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

	return boundaries;
}
} /* namespace natrium */
