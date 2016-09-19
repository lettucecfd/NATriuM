/*
 * DoFBoundary.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_DOFBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_DOFBOUNDARY_H_

namespace natrium {

#include "Boundary.h"

template<size_t dim>
class DoFBoundary: public Boundary<dim> {
private:
	size_t m_boundaryIndicator;

public:
	// TODO boundary indicator
	DoFBoundary() :
			m_boundaryIndicator(0) {
	};

	virtual ~DoFBoundary() {
	};

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	};
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_DOFBOUNDARY_H_ */
