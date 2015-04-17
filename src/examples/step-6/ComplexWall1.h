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

using dealii::Triangulation;

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

private:

	const double m_bottomVelocity;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > makeGrid(double L, size_t refinementLevel);

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
	struct UnstructuredGridFunc
	{
		double refineBoundary (const double y) const {
			return 0.5*(1+std::tanh(y-1));
		}
	  double trans(const double x, const double y) const
	  {
		  double h = 0.7;
		  double k = 6;
	    double upper = 1-h*pow(std::sin(k*4*std::atan(1)*x),2);
	    return y*upper;
	  }
	  dealii::Point<2> operator() (const dealii::Point<2> &in) const
	  {
	    return dealii::Point<2> (in(0),
	                     trans(in(0),refineBoundary(in(1))));
	  }
	};

};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
