/*
 * DensityFluctuation2D.h
 *
 *  Created on: Nov 5, 2014
 *      Author: kk
 */

#ifndef DENSITYFLUCTUATION2D_H_
#define DENSITYFLUCTUATION2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"



namespace natrium {

class DensityFluctuation2D: public ProblemDescription<2> {
public:

	DensityFluctuation2D(double viscosity, size_t refinementLevel);

	virtual ~DensityFluctuation2D();

	virtual void applyInitialDensities(distributed_vector& initialDensities,
			const vector<dealii::Point<2> >& supportPoints) const;

	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<2> >& supportPoints) const;

private:

	shared_ptr<Mesh<2> > makeGrid(size_t refinementLevel);
	shared_ptr<BoundaryCollection<2> > makeBoundaries();

};
} /* namespace natrium */
#endif /* DENSITYFLUCTUATION2D_H_ */
