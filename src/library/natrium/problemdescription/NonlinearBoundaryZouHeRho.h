/*
 * NonlinearBoundaryZouHeRho.h
 *
 *  Created on: 15.12.2015
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_

#include "NonlinearBoundary.h"
#include "BoundaryTools.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

template<size_t dim>
class NonlinearBoundaryZouHeRho {
private:
	boost::shared_ptr<dealii::Function<dim> > m_boundaryPressure;
	size_t m_direction;
	dealii::Tensor<1,dim> m_outwardNormal;
	double m_sign;
public:
	/**
	 * @short
	 * @param[in] direction direction of the outward normal at the wall
	 * 	-# in 2D
	 * 		- 0) left
	 * 		- 1) right
	 * 		- 2) bottom
	 * 		- 3) top
	 * 	-# in 3D
	 * 		- 0) left
	 * 		- 1) right
	 * 		- 2) bottom
	 * 		- 3) top
	 * 		- 4) front
	 * 		- 5) back
	 */
	NonlinearBoundaryZouHeRho(
			boost::shared_ptr<dealii::Function<dim> > boundary_pressure,
			size_t direction);
	virtual ~NonlinearBoundaryZouHeRho() {

	}
	virtual void updateNonlinearBoundaryValues() const;
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_ */
