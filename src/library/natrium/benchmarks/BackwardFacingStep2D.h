/*
 * BackwardFacingStep2D.h
 *
 *  Created on: Jan 5, 2016
 *      Author: kraemer
 */

#ifndef BACKWARDFACINGSTEP2D_H_
#define BACKWARDFACINGSTEP2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"

namespace natrium {

class BackwardFacingStep2D: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
		BackwardFacingStep2D* m_flow;
	public:
		InitialVelocity(BackwardFacingStep2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x,
				const unsigned int component = 0) const;
	};

	class InflowVelocity: public dealii::Function<2> {
	private:
		double m_hStep;
		double m_hDomain;
		double m_uAverage;
	public:
		InflowVelocity(double h_step, double h_domain, double u_av) :
				m_hStep(h_step), m_hDomain(h_domain), m_uAverage(u_av) {
		}
		virtual double value(const dealii::Point<2>& x,
				const unsigned int component = 0) const;
	};

	/// constructor
	BackwardFacingStep2D(double viscosity, double inflow_velocity,
			size_t refinement_level, double L_domain = 18.5,
			double L_step = 2.5, double h_domain = 1.0, double h_step = 0.5);

	/// destructor
	virtual ~BackwardFacingStep2D();

	virtual double getCharacteristicVelocity() const {
		return m_inflowVelocity;
	}

	virtual void refine(Mesh<2>& mesh);

	virtual void transform(Mesh<2>& mesh);
	virtual bool isCartesian(){
		return true;
	}
private:

	const double m_inflowVelocity;
	const double m_LDomain;
	const double m_LStep;
	const double m_HDomain;
	const double m_HStep;
	const size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid();

	/**
	 * @short create boundaries for Couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries(
			double inflow_velocity);

	/**
	 * @short function to generate the unstructured mesh grid
	 */
	/*struct UnstructuredGridFunc {
	 double m_height;
	 double m_ampl;
	 double m_wavelength;
	 double refineBoundaryY(const double y) const {
	 // a smooth curve with gradients that are not 0 or infinity
	 // less steepness on the boundaries leads to fine resolution
	 // of the boundary
	 return 0.5 + 0.5*tanh((y-0.5)*4*atan(1))/tanh((0.5)*4*atan(1));
	 }
	 double refineBoundaryX(const double x) const {
	 // found this formula by analysis with gnuplot;
	 // aim: reduce aspect ratio of cells
	 return  x - 0.019 * m_ampl / 0.49 * sin( 16 * atan(1) * x );
	 }
	 double trans(const double x, const double y) const {
	 double upper = m_height +  m_ampl * std::sin(8 * std::atan(1) * x );
	 return y * upper;
	 }
	 dealii::Point<2> operator()(const dealii::Point<2> &in) const {
	 return dealii::Point<2>(m_wavelength * refineBoundaryX(in(0)), trans(refineBoundaryX(in(0)), refineBoundaryY(in(1))));
	 }
	 UnstructuredGridFunc(double averageHeight, double amplitude, double wavelength) {
	 assert (averageHeight > amplitude);
	 m_height = averageHeight;
	 m_ampl = amplitude;
	 m_wavelength = wavelength;
	 }
	 };*/

};
} /* namespace natrium */

#endif /* BACKWARDFACINGSTEP2D_H_ */
