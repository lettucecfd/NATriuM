/**
 * @file SodShockTube.h
 * @short Description of a simple Periodic Flow (in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SodShockTube_H_
#define SodShockTube_H_

#include "deal.II/grid/tria.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/base/function.h"


#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class SodShockTube: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
		SodShockTube* m_flow;
	public:
		InitialVelocity(SodShockTube* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	class InitialDensity: public dealii::Function<2> {
	private:
		SodShockTube* m_flow;
	public:
		InitialDensity(SodShockTube* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;


	};

	class InitialTemperature: public dealii::Function<2> {
	private:
		SodShockTube* m_flow;
	public:
		InitialTemperature(SodShockTube* flow) :
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
	SodShockTube(int length, double viscosity,
			size_t refinement_level, double u0, double kappa, double perturbation=0.05, double trafo_x=0, double trafo_y=0);

	/// destructor
	virtual ~SodShockTube();

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

			return 1.0/(4.0*M_PI)*(trafo * sin( 4*M_PI*y) + 4*M_PI*y);
		}
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>(trans(in(0),m_tX), trans(in(1),m_tY));
		}
	};


    double m_length;
	double m_u0;
	double m_kappa;
	size_t m_refinementLevel;
	double m_perturbation;
	double m_trafoX;
	double m_trafoY;





	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(int length);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();



};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
