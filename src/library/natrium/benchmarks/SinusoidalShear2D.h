/*
 * SinusoidalShear2D.h
 *
 *  Created on: Feb 2, 2015
 *      Author: kraemer
 */

#ifndef SINUSOIDALSHEAR2D_H_
#define SINUSOIDALSHEAR2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

class SinusoidalShear2D: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
		SinusoidalShear2D* m_flow;
	public:
		InitialVelocity(SinusoidalShear2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};


	/// constructor
	SinusoidalShear2D(double viscosity, double bottomVelocity,
			size_t refinementLevel, double L = 1.0, double averageHeight = 1.0,
			double amplitude = 0.5, double cell_aspect_ratio = 1.0);

	/// destructor
	virtual ~SinusoidalShear2D();

	virtual double getCharacteristicVelocity() const {
		return m_bottomVelocity;
	}

private:

	const double m_bottomVelocity;
	double m_height;
	double m_ampl;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(double cell_aspect_ratio);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries(double bottomVelocity);

	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc {
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
	};

};
} /* namespace natrium */

#endif /* SINUSOIDALGAP2D_H_ */
