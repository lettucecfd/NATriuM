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

#include "../problemdescription/Benchmark.h"
#include "../utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
class TaylorGreenVortex3D: public Benchmark<3> {
public:

  /// constructor
  TaylorGreenVortex3D(double viscosity,
      size_t refinementLevel);

  /// destructor
  virtual ~TaylorGreenVortex3D();

  /**
   * @short analytic solution of the Taylor-Green vortex
   */
  virtual void getAnalyticVelocity(const dealii::Point<3>& x, double t, dealii::Point<3>& velocity) const;


private:

  /**
   * @short create triangulation for couette flow
   * @return shared pointer to a triangulation instance
   */
  shared_ptr<Triangulation<3> > makeGrid(size_t refinementLevel);

  /**
   * @short create boundaries for couette flow
   * @return shared pointer to a vector of boundaries
   * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
   */
  shared_ptr<BoundaryCollection<3> > makeBoundaries();

};

} /* namespace natrium */

#endif /* TAYLORGREENVORTEX3D_H_ */
