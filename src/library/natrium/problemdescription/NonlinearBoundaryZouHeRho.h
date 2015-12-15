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

template <size_t dim>
class NonlinearBoundaryZouHeRho {
private:
	boost::shared_ptr<dealii::Function<dim> > m_boundaryDensity;
public:
	NonlinearBoundaryZouHeRho(boost::shared_ptr<dealii::Function<dim> > boundary_density);
	virtual ~NonlinearBoundaryZouHeRho(){

	}
	virtual void updateNonlinearBoundaryValues() const;
};

} /* namespace natrium */


#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_NONLINEARBOUNDARYZOUHERHO_H_ */
