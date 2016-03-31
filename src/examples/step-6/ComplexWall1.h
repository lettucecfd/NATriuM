/**
 * @file CouetteFlow2D.h
 * @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COMPLEXWALL1_H_
#define COMPLEXWALL1_H_

#include "deal.II/grid/tria.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Couette Flow (regular channel flow in square domain).
 *  The domain is [0,1]^2. The top plate is moved with constant velocity. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 *  @note The analytic solution is obtained by a formula stated in Min and Lee (2011).
 */
class ComplexWall1: public ProblemDescription<2> {
public:

	/// constructor
	ComplexWall1(double viscosity, double bottomVelocity , size_t refinementLevel, double L=1.0);

	/// destructor
	virtual ~ComplexWall1();

	virtual double getCharacteristicVelocity() const {
		return m_bottomVelocity;
	}

	virtual void refine();

	virtual void transform(Mesh<2>& mesh);

private:

	const double m_bottomVelocity;
	const size_t m_refinementLevel;
	const double m_L;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(double L, size_t refinementLevel);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries(double bottomVelocity);

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
	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc {
		double m_L;
		double refineBoundaryY(const double y) const {
			// a smooth curve with gradients that are not 0 or infinity
			// less steepness on the boundaries leads to fine resolution
			// of the boundary
			return 1 * (0.5 + 0.5*tanh((y / 1 -0.5)*4*atan(1))/tanh((0.5)*4*atan(1)));
		}
		double refineBoundaryX(const double x) const {
			// found this formula by analysis with gnuplot;
			// aim: reduce aspect ratio of cells
			//return  x - 0.019 * m_ampl / 0.49 * sin( 16 * atan(1) * x );
			// dx ~ local height
			return x;
		}
		double trans(const double x, const double y) const {
			double upper = 1 +  0.3 * cos(8 * atan(1) * x / (m_L*1)) + 0.1 * pow(cos(8.0 * atan(1) * x / (m_L / 6.0)),3) + 0.4 * cos(1 + 8.0 * atan(1) * x / (m_L / 2.0)) + 0.1 * cos(0.12 + 8.0 * atan(1) * x / (m_L / 5.0));

			return y * upper;
		}
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>(refineBoundaryX(in(0)), trans(refineBoundaryX(in(0)), refineBoundaryY(in(1))));
		}
		UnstructuredGridFunc(double L){
			m_L = L;
		}
	};

};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
