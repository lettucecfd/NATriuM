/*
 * DoFBoundary.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_DOFBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_DOFBOUNDARY_H_

#include "Boundary.h"

#include "deal.II/base/function.h"

#include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim>
class DoFBoundary: public Boundary<dim> {
private:
	size_t m_boundaryIndicator;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryDensity;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryVelocity;

public:

	DoFBoundary(size_t boundaryIndicator) :
			m_boundaryIndicator(boundaryIndicator) {

	}
	// TODO boundary indicator
	DoFBoundary(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
			m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
					boundaryDensity), m_boundaryVelocity(boundaryVelocity) {
	}
	;

	virtual ~DoFBoundary() {
	}
	;

	virtual bool isPeriodic() const {
		return false;
	}
	virtual bool isLinearFluxBoundary() const {
		return false;
	}
	virtual bool isDoFBoundary() const {
		return true;
	}

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}
	const boost::shared_ptr<dealii::Function<dim> >& getBoundaryDensity() const {
		return m_boundaryDensity;
	}

	const boost::shared_ptr<dealii::Function<dim> >& getBoundaryVelocity() const {
		return m_boundaryVelocity;
	}

	;
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_DOFBOUNDARY_H_ */
