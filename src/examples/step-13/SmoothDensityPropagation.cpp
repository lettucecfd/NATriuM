/**
 * @file SmoothDensityPropagation.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SmoothDensityPropagation.h"

#include "deal.II/grid/grid_generator.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "natrium/boundaries/PeriodicBoundary.h"

#include "natrium/utilities/Math.h"

namespace natrium {

    SmoothDensityPropagation::SmoothDensityPropagation(double viscosity,
                                                       size_t refinementLevel, double cs, bool init_rho_analytically,
                                                       double L , double horizontal) :
            CompressibleBenchmark<2>(makeGrid(L), viscosity, L), m_cs(
            cs), m_analyticInit(init_rho_analytically), m_refinementLevel(refinementLevel), m_horizontalVelocity(horizontal),
            m_L(L) {

        /// apply boundary values
        setBoundaries(makeBoundaries());
        // apply initial and analytical solution
        this->setAnalyticU(boost::make_shared<AnalyticVelocity>(this));
        this->setAnalyticRho(boost::make_shared<AnalyticDensity>(this));
        this->setAnalyticT(boost::make_shared<AnalyticTemperature>(this));
        assert (L > 0);
    }

    SmoothDensityPropagation::~SmoothDensityPropagation() {
    }

    double SmoothDensityPropagation::AnalyticTemperature::value(const dealii::Point<2> &x,
                                                                const unsigned int component) const {
        assert(component < 2);
        return 1.0 / this->m_flow->getAnalyticRhoFunction(this->get_time())->value(x,component);
    }

    double SmoothDensityPropagation::AnalyticVelocity::value(const dealii::Point<2> &x,
                                                             const unsigned int component) const {
        assert(component < 2);
        if (component == 0) {
            return m_flow->getHorizontalVelocity();
        } else {
            return 0.0;
        }
    }

    double SmoothDensityPropagation::AnalyticDensity::value(const dealii::Point<2> &x,
                                                            const unsigned int component) const {
        assert (component == 0);

        double rho0 = 1.0;
        double rho = rho0 + 0.2*(sin(2.0 * M_PI * x(0) / m_flow->m_L - 2*M_PI*m_flow->getHorizontalVelocity() * this->get_time()) *
                            sin(x(1) * (2 * M_PI) / m_flow->m_L));
        return rho;

    }

/**
 * @short create triangulation for couette flow
 * @return shared pointer to a triangulation instance
 */
    boost::shared_ptr<Mesh<2> > SmoothDensityPropagation::makeGrid(double L) {
        //Creation of the principal domain
#ifdef WITH_TRILINOS_MPI
        boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >(MPI_COMM_WORLD);
#else
        boost::shared_ptr<Mesh<2> > square = boost::make_shared<Mesh<2> >();
#endif
        dealii::GridGenerator::hyper_cube(*square, 0, L);// * atan(1));

        // Assign boundary indicators to the faces of the "parent cell"
        Mesh<2>::active_cell_iterator cell = square->begin_active();
        cell->face(0)->set_all_boundary_ids(0);  // left
        cell->face(1)->set_all_boundary_ids(1);  // right
        cell->face(2)->set_all_boundary_ids(2);  // top
        cell->face(3)->set_all_boundary_ids(3);  // bottom


        return square;
    }

/**
 * @short create boundaries for couette flow
 * @return shared pointer to a vector of boundaries
 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
 */
    boost::shared_ptr<BoundaryCollection<2> > SmoothDensityPropagation::makeBoundaries() {

        // make boundary description
        boost::shared_ptr<BoundaryCollection<2> > boundaries = boost::make_shared<
                BoundaryCollection<2> >();
        boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(0, 1, 0, getMesh()));
        boundaries->addBoundary(boost::make_shared<PeriodicBoundary<2> >(2, 3, 1, getMesh()));

        // Get the triangulation object (which belongs to the parent class).
        boost::shared_ptr<Mesh<2> > tria_pointer = getMesh();

        return boundaries;
    }
} /* namespace natrium */
