/*
 * TaylorGreenVortex3D.h
 *
 *  Created on: Sep 18, 2014
 *      Author: bajat
 */

#ifndef TAYLORGREENVORTEX3D_H_
#define TAYLORGREENVORTEX3D_H_


/**
 * @file TaylorGreenVortex3D.h
 * @short Description of a simple Periodic Flow (in cubic domain).
 */

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class TaylorGreenVortex3D: public ProblemDescription<3> {
public:

	/**
	 * @short class to describe the x-component of the initial velocity
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class InitialVelocity: public dealii::Function<3> {
	private:
		TaylorGreenVortex3D* m_flow;
		double m_U;
	public:
		InitialVelocity(TaylorGreenVortex3D* flow, double u) :
				m_flow(flow), m_U(u) {
		}
		virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;
	};

	/**
	 * @short class to describe the initial density (i.e. the initial pressure in LBM)
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class InitialDensity: public dealii::Function<3> {
	private:
		TaylorGreenVortex3D* m_flow;
		double m_csSquare;
	public:
		InitialDensity(TaylorGreenVortex3D* flow, double cs) :
				m_flow(flow), m_csSquare(cs*cs) {
		}
		virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;
	};

  /// constructor
	TaylorGreenVortex3D(double viscosity,
			size_t refinementLevel, double u, double cs);

  /// destructor
  virtual ~TaylorGreenVortex3D();

private:

  /**
   * @short create triangulation for couette flow
   * @return shared pointer to a triangulation instance
   */
  boost::shared_ptr<Mesh<3> > makeGrid();

  /**
   * @short create boundaries for couette flow
   * @return shared pointer to a vector of boundaries
   * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
   */
  boost::shared_ptr<BoundaryCollection<3> > makeBoundaries();

};

} /* namespace natrium */

#endif /* TAYLORGREENVORTEX3D_H_ */
