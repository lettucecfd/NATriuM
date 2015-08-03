/*
 * LubrcationSine.h
 *
 *  Created on: Feb 2, 2015
 *      Author: kraemer
 */

#ifndef LUBRICATIONSINE_H_
#define LUBRICATIONSINE_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

class LubricationSine: public ProblemDescription<2> {
public:

	/// constructor
	LubricationSine(double viscosity, double bottomVelocity,
			size_t refinementLevel, double L = 1.0, double averageHeight = 1.0,
			double amplitude = 0.5, double cellAspectRatio=0.5, double roughnessHeight = 0.0, size_t roughnessLengthRatio = 1);

	/// destructor
	virtual ~LubricationSine();

	virtual double getCharacteristicVelocity() const {
		return m_bottomVelocity;
	}

private:

	const double m_bottomVelocity;
	double m_height;
	double m_ampl;
	const double m_length;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Mesh<2> > makeGrid(double L, size_t refinementLevel,
			double averageHeight, double amplitude, double cellAspectRatio,  double roughnessHeight, size_t roughnessLengthRatio);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	shared_ptr<BoundaryCollection<2> > makeBoundaries(double bottomVelocity);

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialDensities(distributed_vector& initialDensities,
			const vector<dealii::Point<2> >& supportPoints) const;

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<2> >& supportPoints) const;

	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc {
		double m_height;
		double m_ampl;
		double m_wavelength;
		double m_roughnessHeight;
		size_t m_roughnessLengthRatio;
		double refineBoundaryY(const double y) const {
			// a smooth curve with gradients that are not 0 or infinity
			// less steepness on the boundaries leads to fine resolution
			// of the boundary
			return m_height * (0.5 + 0.5*tanh((y / m_height -0.5)*4*atan(1))/tanh((0.5)*4*atan(1)));
		}
		double refineBoundaryX(const double x) const {
			// found this formula by analysis with gnuplot;
			// aim: reduce aspect ratio of cells
			//return  x - 0.019 * m_ampl / 0.49 * sin( 16 * atan(1) * x );
			// dx ~ local height
			return x + m_wavelength * m_ampl / m_height / (8.0 * atan(1.0)) * sin (8.0 * atan(1.0) * x / m_wavelength);
		}
		double trans(const double x, const double y) const {
			double upper = m_height +  m_ampl * cos(8 * atan(1) * x / m_wavelength) + m_roughnessHeight * cos(8.0 * atan(1) * x / (m_wavelength/m_roughnessLengthRatio));
			return y / m_height * upper;
		}
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>(refineBoundaryX(in(0)), trans(refineBoundaryX(in(0)), refineBoundaryY(in(1))));
		}
		UnstructuredGridFunc(double averageHeight, double amplitude, double wavelength, double roughnessHeight, size_t roughnessLengthRatio) {
			assert (averageHeight > amplitude);
			m_height = averageHeight;
			m_ampl = amplitude;
			m_wavelength = wavelength;
			m_roughnessHeight = roughnessHeight;
			m_roughnessLengthRatio = roughnessLengthRatio;
		}
	};

};
} /* namespace natrium */

#endif /* SINUSOIDALGAP2D_H_ */
