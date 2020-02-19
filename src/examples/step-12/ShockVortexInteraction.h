/**
 * @file ShockVortexInteraction.h
 * @short Shock-vortex interaction according to Inoue (1999)
 * @date 14.08.2019
 * @author Dominik Wilde, University of Siegen
 */

#ifndef ShockVortexInteraction_H_
#define ShockVortexInteraction_H_

#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/base/function.h"


#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"
#include <tuple>



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class ShockVortexInteraction: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
                ShockVortexInteraction* m_flow;
	public:
                InitialVelocity(ShockVortexInteraction* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	class InitialDensity: public dealii::Function<2> {
	private:
                ShockVortexInteraction* m_flow;
	public:
                InitialDensity(ShockVortexInteraction* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;


	};

	class InitialTemperature: public dealii::Function<2> {
	private:
                ShockVortexInteraction* m_flow;
	public:
                InitialTemperature(ShockVortexInteraction* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};


	/// constructor (with default cs=1/sqrt(3))
	/**
	 * @short Constructor 
	 * @param trafo_x < 1: 0=regular grid spacing
	 * @param trafo_y < 1: 0=regular grid spacing
	 **/
        ShockVortexInteraction(double viscosity,
			size_t refinement_level, double u0, double kappa, double Ma_v, double perturbation=0.05, double trafo_x=0, double trafo_y=0);

	/// destructor
        virtual ~ShockVortexInteraction();

	virtual void refine(Mesh<2>& mesh){
		// Refine grid
		mesh.refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<2>& mesh){
		dealii::GridTools::transform(
						UnstructuredGridFunc(m_trafoX, m_trafoY), mesh);
		std::ofstream out_file("/tmp/grid_out.eps");
		dealii::GridOut().write_eps(*getMesh(), out_file);
		out_file.close();
	}
	virtual bool isCartesian(){
		return true;
	}

    std::tuple<double,double> calcOffsets(double left, double right){

        double d = (right/2.0-left/2.0);
        double c = 2.0*(left)/(right-left)+1.0;
        return {c,d};
    }
private:

	struct UnstructuredGridFunc {
		double m_tX;
		double m_tY;
		UnstructuredGridFunc(double trafo_x, double trafo_y):
	       	m_tX (trafo_x), m_tY (trafo_y)	{
		}
		double trans(const double y, double trafo) const {

            return 36./(M_PI)*(trafo * sin(M_PI/36.*(y+6.0)) + M_PI/36.*y);
		}
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>(trans(in(0),m_tX), trans(in(1),m_tY));
		}
	};


	double m_u0;
	double m_kappa;
	size_t m_refinementLevel;
	double m_perturbation;
	double m_trafoX;
	double m_trafoY;

	double m_Ma_v;
	double m_shockPosition = 30;
	double m_shockSteepness = 1000000;
    double m_densityLeft = 1.34161490;
    double m_densityRight = 1.0;
    double m_machNumberLeft = 0.84217047;
	double m_machNumberRight = 1.2;
    double m_temperatureLeft = 1.12799382716;
    double m_temperatureRight = 1.0;




	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid();

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();



};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
