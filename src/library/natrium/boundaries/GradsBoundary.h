/*
 * GradsBoundaryCondition.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_

#include <array>
#include "../stencils/Stencil.h"

namespace natrium {

template<size_t dim>
class GradsBoundary: public DoFBoundary<dim> {
public:
	GradsBoundary() {

	}
	virtual ~GradsBoundary() {

	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSBOUNDARYCONDITION_H_ */
